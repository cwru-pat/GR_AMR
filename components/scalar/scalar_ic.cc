#include "../../cosmo_includes.h"
#include "scalar_ic.h"

using namespace SAMRAI;

namespace cosmo
{

void scalar_ic_set_semianalytic_test(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  boost::shared_ptr<tbox::Database> cosmo_scalar_db)
{

  boost::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;


  const double * domain_lower = &grid_geometry.getXLower()[0];
  const double * domain_upper = &grid_geometry.getXUpper()[0];

  real_t L[3];
  
  for(int i = 0 ; i < DIM; i++)
    L[i] = domain_upper[i] - domain_lower[i];

  const double * dx = &grid_geometry.getDx()[0];
  
  idx_t NX = round(L[0] / dx[0]);
  //idx_t NY = round(L[1] / dx[1]);
  //idx_t NZ = round(L[2] / dx[2]);

  fftw_complex *f_temp;  
  f_temp = (fftw_complex *) fftw_malloc((NX/2+1)
                                         *((long long) sizeof(fftw_complex)));

  
  double *r_temp = new double[NX], *der_bak = new double[NX];
  memset(r_temp, 0, NX*sizeof(double));


  real_t Lambda = 1.0;
  
  for(int i = 0; i < NX; i++)
  {
    real_t phi_temp = 1.0 + 0.01 * std::sin(4.0 * PI *( (real_t)i / NX - 0.125));
    real_t K_temp = 9.763423957014197;
    if(2.0 * i <= NX)
      r_temp[i] = std::sqrt(std::fabs(
                              (pw2(4.0 * PI) * 0.01  *  std::sin(4.0 * PI *( (real_t)i / NX - 0.125) )
                               +(-2.0 * PI * Lambda + pw2(K_temp)/12.0 )*
                               std::pow(phi_temp, 5.0) )
                              /(PI * phi_temp ))) ;
    else
      r_temp[i]= -std::sqrt(std::fabs(
                              (pw2(4.0 * PI) * 0.01  *  std::sin(4.0 * PI *( (real_t)i / NX - 0.125) )
                               +(-2.0 * PI * Lambda + pw2(K_temp)/12.0 )*
                               std::pow(phi_temp, 5.0) )
                              /(PI * phi_temp ))) ;
    der_bak[i] = r_temp[i];
  }

  // plans for taking FFTs
  fftw_plan p_c2r;
  fftw_plan p_r2c;

  p_c2r = fftw_plan_dft_c2r_1d(NX,
                               f_temp, r_temp,
                               FFTW_ESTIMATE);
  p_r2c = fftw_plan_dft_r2c_1d(NX,
                               r_temp, f_temp,
                               FFTW_ESTIMATE);

  fftw_execute(p_r2c);

  for(int i = 0; i < NX/2 + 1; i++)
  {
    if(i > 0)
    {
      std::swap(f_temp[i][0], f_temp[i][1]);
      f_temp[i][0] = dx[0] * f_temp[i][0] /
        ( 2.0 * PI * (real_t) i / NX);
      f_temp[i][1] = - dx[0] * f_temp[i][1] /
        ( 2.0 * PI * (real_t) i / NX);

    }
    else
      f_temp[i][0] = f_temp[i][1] = 0;
  }

  fftw_execute(p_c2r);
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;

    bssn->initPData(patch);
    bssn->initMDA(patch);

    scalar->initPData(patch);
    scalar->initMDA(patch);
    
    arr_t & DIFFchi = bssn->DIFFchi_a;
    arr_t & phi = scalar->phi_a; // field
    arr_t & psi1 = scalar->psi1_a; // derivative of phi in x-dir
    arr_t & psi2 = scalar->psi2_a; // derivative of phi in y-dir
    arr_t & psi3 = scalar->psi3_a; // derivative of phi in z-dir
  
    arr_t & K = bssn->DIFFK_a; // extrinsic curvature


    const hier::Box& box = bssn->DIFFchi_a_pdata->getGhostBox();
    const hier::Box& inner_box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    const int * inner_lower = &inner_box.lower()[0];
    const int * inner_upper = &inner_box.upper()[0];

    
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          idx_t ii = (i + NX ) % NX; 
          phi(i, j, k) = r_temp[ii] / NX;
          DIFFchi(i, j, k) = 1.0 / pw2(1.0 + 0.01 * std::sin(4.0 * PI *( (real_t)ii / NX - 0.125))) - 1.0;
          K(i, j, k) = 9.763423957014197;
        }
      }
    }


    for(int k = inner_lower[2]; k <= inner_upper[2]; k++)
    {
      for(int j = inner_lower[1]; j <= inner_upper[1]; j++)
      {
        for(int i = inner_lower[0]; i <= inner_upper[0]; i++)
        {
          psi1(i, j, k) = derivative(i, j, k, 1, phi, dx);
          psi2(i, j, k) = derivative(i, j, k, 2, phi, dx);
          psi3(i, j, k) = derivative(i, j, k, 3, phi, dx);
        }
      }
    }   
  }
  
}

}
