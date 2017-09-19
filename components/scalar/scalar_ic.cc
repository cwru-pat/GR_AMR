#include "../../cosmo_includes.h"
#include "scalar_ic.h"
#include "SAMRAI/tbox/HDFDatabase.h"
#include "../elliptic_solver/full_multigrid.h"
#include "../elliptic_solver/multigrid_bd_handler.h"
#include "../../utils/Array.h"
#include <iostream>
#include <fstream>

#include "../horizon/AHFD/jtutil/util.hh"
#include "../horizon/AHFD/jtutil/array.hh"
#include "../horizon/AHFD/jtutil/util_String.h"
#include "../horizon/AHFD/jtutil/util_Table.h"
#include "../horizon/AHFD/jtutil/interpolator/InterpLocalUniform.h"
//#include "../horizon/AHFD/AHFD_types.h"
#include "../horizon/AHFD/AHFD_macros.h"

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

  /*
    test for a spherical gaussian field collapse under 
    sommerfield boundary condition
  */
bool scalar_ic_set_scalar_collapse_sommerfield(
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
  double dx[3];
  
  for(int i = 0 ; i < DIM; i++)
  {
    L[i] = domain_upper[i] - domain_lower[i];
    dx[i] = (grid_geometry.getDx()[i]) / (1<<ln);
  }

  double rho_tot = 0;
  std::ifstream myReadFile;
  myReadFile.open("init.dat");

  double temp2[50001];

  int cnt =0;
  if (myReadFile.is_open()) {
    while (!myReadFile.eof()) {


      myReadFile >> temp2[cnt];
      cnt++;

    }
  }
    
  double r0 = cosmo_scalar_db->getDoubleWithDefault("r0", 0);
  double sigma = cosmo_scalar_db->getDoubleWithDefault("sigma", 1.0);
  double q = cosmo_scalar_db->getDoubleWithDefault("q", 2.0);
  real_t phi_0 = cosmo_scalar_db->getDoubleWithDefault("phi_0", 1.0);

  int n_max = cosmo_scalar_db->getIntegerWithDefault("n_max", 1);

  real_t delta_phi = cosmo_scalar_db->getDoubleWithDefault("delta_phi", 0.1);


  /********/

  //  initializing interpolator

 /**********/

  char interpolator_par[] = "order=3 boundary_off_centering_tolerance={1.0e-10 1.0e-10} boundary_extrapolation_tolerance={0.0 0.0}";

  int interp_handle =  Util_TableCreateFromString(interpolator_par);

  if(interp_handle < 0)
    TBOX_ERROR("Failed to creat interp table!\n");
  
  const double origin[1] = {0.0005};
  const double delta[1] = {0.001};

  const void* input_arrays[1];

  const CCTK_INT input_array_dims[1] = {50000};
  const CCTK_INT input_array_type_codes[1] = {CCTK_VARIABLE_REAL};

  CCTK_REAL    output_real[1];
  const CCTK_INT output_array_type_codes[1] = {CCTK_VARIABLE_REAL};  

  const void* interp_coords[1];
  void* output_arrays[1];
  
  input_arrays [0] = (const void *) temp2;  
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;

    bssn->initPData(patch);
    bssn->initMDA(patch);

    scalar->initPData(patch);
    scalar->initMDA(patch);
    
    arr_t & DIFFchi_a = bssn->DIFFchi_a;
    arr_t & phi_a = scalar->phi_a; // field
    arr_t & psi1_a = scalar->psi1_a; // derivative of phi in x-dir
    arr_t & psi2_a = scalar->psi2_a; // derivative of phi in y-dir
    arr_t & psi3_a = scalar->psi3_a; // derivative of phi in z-dir
  
    arr_t & K_a = bssn->DIFFK_a; // extrinsic curvature

    arr_t & DIFFalpha_a = bssn->DIFFalpha_a;
    arr_t & DIFFgamma11_a = bssn->DIFFgamma11_a;
    arr_t & DIFFgamma22_a = bssn->DIFFgamma22_a;
    arr_t & DIFFgamma33_a = bssn->DIFFgamma33_a;
    arr_t & DIFFgamma12_a = bssn->DIFFgamma12_a;
    arr_t & DIFFgamma13_a = bssn->DIFFgamma13_a;
    arr_t & DIFFgamma23_a = bssn->DIFFgamma23_a;
    
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
          K_a(i, j, k) = 0;

          double x = ((double)i + 0.5) * dx[0] - L[0] / 2.0;
          double y = ((double)j + 0.5) * dx[1] - L[1] / 2.0;
          double z = ((double)k + 0.5) * dx[2] - L[2] / 2.0;

          double r = sqrt(x * x + y * y + z * z);

          double st = sqrt(x * x + y * y) / r;

          const CCTK_REAL interp_r[1] = {r};

          CCTK_REAL output[1];
          
          interp_coords[0] = (const void *) interp_r;
          output_arrays[0] = (      void *) output;
          if (AEILocalInterp_U_Hermite(1,  
                                       interp_handle,  
                                       origin, delta,  
                                       1,  
                                       CCTK_VARIABLE_REAL,  
                                       interp_coords,  
                                       1,  
                                       input_array_dims,  
                                       input_array_type_codes,  
                                       input_arrays,  
                                       1,  
                                       output_array_type_codes,  
                                       output_arrays) < 0)
            TBOX_ERROR("error return from interpolator!");

          
          //          double aa = (temp2[l] * ((double)(l+2) * 0.001 - r) + temp2[l+1] * (r - (double)(l+1)*0.001) ) / 0.001;
          double aa = output[0];


          DIFFgamma11_a(i, j, k) = (pw2(aa)*pw2(x) + pw2(y) + pw2(z))/(pw2(x) + pw2(y) + pw2(z)) - 1.0;

          DIFFgamma22_a(i, j, k) = (pw2(x) + pw2(aa)*pw2(y) + pw2(z))/(pw2(x) + pw2(y) + pw2(z)) - 1.0;

          DIFFgamma33_a(i, j, k) = (pw2(x) + pw2(y) + pw2(aa)*pw2(z))/(pw2(x) + pw2(y) + pw2(z)) - 1.0;

          DIFFgamma12_a(i, j, k) = ((-1 + pw2(aa))*x*y)/(pw2(x) + pw2(y) + pw2(z));
          DIFFgamma13_a(i, j, k) = ((-1 + pw2(aa))*x*z)/(pw2(x) + pw2(y) + pw2(z));

          DIFFgamma23_a(i, j, k) = ((-1 + pw2(aa))*y*z)/(pw2(x) + pw2(y) + pw2(z));

          DIFFchi_a(i, j, k) = 0;
          phi_a(i, j, k) = phi_0 + delta_phi * /*(1.0+pw2(r)) */ pw2(r) *
            exp( - pow(fabs( (r - r0) / sigma) , q)) ;
            
          DIFFalpha_a(i, j, k) = 0;
        }
      }
    }

    for(int k = inner_lower[2]; k <= inner_upper[2]; k++)
    {
      for(int j = inner_lower[1]; j <= inner_upper[1]; j++)
      {
        for(int i = inner_lower[0]; i <= inner_upper[0]; i++)
        {
          psi1_a(i, j, k) = derivative(i, j, k, 1, phi_a, dx);
          psi2_a(i, j, k) = derivative(i, j, k, 2, phi_a, dx);
          psi3_a(i, j, k) = derivative(i, j, k, 3, phi_a, dx);              
        }
      }
    }   
  }

  Util_TableDestroy(interp_handle);
  
  return true;


}
  
bool scalar_ic_set_scalar_collapse(
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
  double dx[3];
  
  for(int i = 0 ; i < DIM; i++)
  {
    L[i] = domain_upper[i] - domain_lower[i];
    dx[i] = (grid_geometry.getDx()[i]) / (1<<ln);
  }

  std::string boundary_type = "periodic";
  multigridBdHandler * bd_handler = new multigridBdHandler(boundary_type, L, 10);
  
  idx_t NX = round(L[0] / dx[0]);
  idx_t NY = round(L[1] / dx[1]); 
  idx_t NZ = round(L[2] / dx[2]);

  
  /******getting some parameters from input database****/

  real_t phi_0 = cosmo_scalar_db->getDoubleWithDefault("phi_0", 1.0);

  int n_max = cosmo_scalar_db->getIntegerWithDefault("n_max", 1);

  real_t delta_phi = cosmo_scalar_db->getDoubleWithDefault("delta_phi", 0.1);

  // solve for BSSN fields using multigrid class:
  real_t relaxation_tolerance = cosmo_scalar_db->getDoubleWithDefault("relaxation_tolerance", 1e-8);

  int num_vcycles = cosmo_scalar_db->getIntegerWithDefault("vcycles", 20);

  double DIFFalpha_0 = cosmo_scalar_db->getDoubleWithDefault("DIFFalpha", 0);

  /******ending collecting parameters from input database****/

  //  double * phi = new double[(NX+2*STENCIL_ORDER) * (NY+2*STENCIL_ORDER) * (NZ+2*STENCIL_ORDER)];
  CosmoArray<idx_t, real_t>  phi;
  phi.init(NX, NY, NZ);
  //  std::vector<double> phi(NX*NY*NZ, phi_0);
  
  LOOP3()
    phi[INDEX(i, j, k)] = phi_0;

  std::random_device rd;
  std::mt19937 gen(7.0 /*rd()*/);
  std::uniform_real_distribution<real_t> dist(0, 2.0*PI);


  std::string initial_type =
    cosmo_scalar_db->getStringWithDefault("initial_type", "modes");
  if(initial_type == "modes")
  {
    for(int n = -n_max; n <= n_max; ++n)
    {
      if(n != 0)
      {
        // random phases
        real_t x_phase = dist(gen),
          y_phase = dist(gen),
          z_phase = dist(gen);

        LOOP3()
        {
          // some sinusoidal modes
          phi[INDEX(i,j,k)] += delta_phi*(
            cos(2.0*PI*((real_t) n/NX)*((real_t)i + 0.5) )
            + cos(2.0*PI*((real_t) n/NY)*((real_t)j + 0.5) )
            + cos(2.0*PI*((real_t) n/NZ)*((real_t)k+ 0.5) )
          );
        }
      }
    }
  }
  else if(initial_type == "gaussian")
  {
    double r0 = cosmo_scalar_db->getDoubleWithDefault("r0", 0);
    double sigma = cosmo_scalar_db->getDoubleWithDefault("sigma", 1.0);
    double q = cosmo_scalar_db->getDoubleWithDefault("q", 2.0);
    LOOP3()
    {
      double x = L[0] / NX * ((double)i + 0.5) - L[0] / 2.0;
      double y = L[1] / NX * ((double)j + 0.5) - L[1] / 2.0;
      double z = L[2] / NX * ((double)k + 0.5) - L[2] / 2.0;
      double r = sqrt(x * x + y * y + z * z);
      phi[INDEX(i,j,k)] = delta_phi * pw3(r) *
      exp( - pow(fabs( (r - r0) / sigma) , q)) ;
    }
  }
    
  bd_handler->fillBoundary(phi._array, phi.nx, phi.ny, phi.nz);
  // compute background/average K
  real_t K_src = 0;
  LOOP3()
  {
    BSSNData bd = {0};
    ScalarData sd = {0};
    
    sd.phi = phi[INDEX(i,j,k)];
    
    K_src += (pw2((1.0/12.0*phi[INDEX(i-2,j,k)] - 2.0/3.0*phi[INDEX(i-1,j,k)] + 2.0/3.0*phi[INDEX(i+1,j,k)]- 1.0/12.0*phi[INDEX(i+2,j,k)])/dx[0])
              + pw2((1.0/12.0*phi[INDEX(i,j-2,k)] - 2.0/3.0*phi[INDEX(i,j-1,k)] + 2.0/3.0*phi[INDEX(i,j+1,k)]- 1.0/12.0*phi[INDEX(i,j+2,k)])/dx[1])
              +pw2((1.0/12.0*phi[INDEX(i,j,k-2)] - 2.0/3.0*phi[INDEX(i,j,k-1)] + 2.0/3.0*phi[INDEX(i,j,k+1)]- 1.0/12.0*phi[INDEX(i,j,k+2)])/dx[2])) / 2.0
      + scalar->potentialHandler->ev_potential(&bd, &sd);
  }
  
  K_src = -std::sqrt(24.0 * PI * K_src/NX/NY/NZ);
  
  
  boost::shared_ptr<tbox::HDFDatabase > hdf (new tbox::HDFDatabase("hdf_db"));

  std::string filename = "h5_data_lv_";

  filename += tbox::Utilities::intToString(ln, 1);

  CosmoArray<idx_t, real_t> * DIFFchi = new CosmoArray<idx_t, real_t> [1];

  bool flag = false;

  DIFFchi[0].init(NX, NY, NZ);

  if(exist(filename))
  {
    hdf->open(filename, 1);
    const std::vector<double> & temp = hdf->getDoubleVector("DIFFchi");

    // if file exist but corresponding database not exist
    if(temp.empty())
      TBOX_ERROR("Getting empty array from file "<<filename<<"\n");

    tbox::pout<<"Read initial configuration database for level "<<ln<<"\n";
    
    for(int i = 0; i < temp.size(); i++)
      DIFFchi[0]._array[i] = temp[i];
    
    flag = true;
  }
  else
  {
    // create and open the file
    hdf->create(filename);
    hdf->open(filename, 1);

    idx_t molecule_n[] = {3};
    
    FASMultigrid multigrid(
      DIFFchi, 1, molecule_n, 4, 5, relaxation_tolerance, L, NX, NY, NZ, bd_handler);

    atom atom_tmp = {0};

    //initializing equations
    multigrid.eqns[0][0].init(1, 1);
    multigrid.eqns[0][1].init(1, 1);
    multigrid.eqns[0][2].init(1, 1);


    //adding terms to eqn
    //add first laplacian term
    atom_tmp.type = multigrid.atom_type::lap;
    atom_tmp.u_id = 0;
    multigrid.eqns[0][0].add_atom(atom_tmp);

    //add second term
  
    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 1;
    multigrid.eqns[0][1].add_atom(atom_tmp);

    //add third term

    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 5;
    multigrid.eqns[0][2].add_atom(atom_tmp);

    real_t avg1 = 0.0, avg5 = 0.0;
    LOOP3()
    {
      BSSNData bd = {0};
      ScalarData sd = {0};
    
      sd.phi = phi[INDEX(i,j,k)];

      real_t value = PI*
        (pw2((1.0/12.0*phi[INDEX(i-2,j,k)] - 2.0/3.0*phi[INDEX(i-1,j,k)] + 2.0/3.0*phi[INDEX(i+1,j,k)]- 1.0/12.0*phi[INDEX(i+2,j,k)])/dx[0])
         + pw2((1.0/12.0*phi[INDEX(i,j-2,k)] - 2.0/3.0*phi[INDEX(i,j-1,k)] + 2.0/3.0*phi[INDEX(i,j+1,k)]- 1.0/12.0*phi[INDEX(i,j+2,k)])/dx[1])
         +pw2((1.0/12.0*phi[INDEX(i,j,k-2)] - 2.0/3.0*phi[INDEX(i,j,k-1)] + 2.0/3.0*phi[INDEX(i,j,k+1)]- 1.0/12.0*phi[INDEX(i,j,k+2)])/dx[2]));

      avg1 += value;
      multigrid.setPolySrcAtPt(0, 1, i, j, k, value); //set value for term 1

      value = 2.0* PI*scalar->potentialHandler->ev_potential(&bd, &sd) - K_src * K_src / 12.0;
      multigrid.setPolySrcAtPt(0, 2, i, j, k, value); //set value for term 2
      avg5 += value;
    }
    avg1 = avg1/NX/NY/NZ;
    avg5 = avg5/NX/NY/NZ;
    multigrid.initializeRhoHeirarchy();

    if(avg1 * avg5 > 0)
      TBOX_ERROR("Cannot find proper initial setting for phi\n");
    
    LOOP3()
    {
      idx_t idx = INDEX(i,j,k);
      double x = ((double)i + 0.5) * dx[0] - L[0] / 2.0;
      double y = ((double)j + 0.5) * dx[1] - L[1] / 2.0;
      double z = ((double)k + 0.5) * dx[2] - L[2] / 2.0;

      double r = sqrt(x * x + y * y + z * z);
      
      DIFFchi[0][idx] = std::pow(-avg1/avg5,1.0/4.0);
    }
    //    std::cout<<std::pow(-avg1/avg5,1.0/4.0)<<"\n";
    bd_handler->fillBoundary(DIFFchi[0]._array, DIFFchi[0].nx, DIFFchi[0].ny, DIFFchi[0].nz);

    multigrid.VCycles(num_vcycles);

    LOOP3()
    {
      idx_t idx = INDEX(i, j, k);
      DIFFchi[0][idx] = 1.0 / pw2(DIFFchi[0][idx]) - 1.0;
    }
    bd_handler->fillBoundary(DIFFchi[0]._array, DIFFchi[0].nx, DIFFchi[0].ny, DIFFchi[0].nz);

    hdf->putDoubleArray("DIFFchi", DIFFchi[0]._array, (NX+2*STENCIL_ORDER)*(NY+2*STENCIL_ORDER)*(NZ+2*STENCIL_ORDER));
    flag = true;
  }

  hdf->close();
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;

    bssn->initPData(patch);
    bssn->initMDA(patch);

    scalar->initPData(patch);
    scalar->initMDA(patch);
    
    arr_t & DIFFchi_a = bssn->DIFFchi_a;
    arr_t & phi_a = scalar->phi_a; // field
    arr_t & psi1_a = scalar->psi1_a; // derivative of phi in x-dir
    arr_t & psi2_a = scalar->psi2_a; // derivative of phi in y-dir
    arr_t & psi3_a = scalar->psi3_a; // derivative of phi in z-dir
  
    arr_t & K_a = bssn->DIFFK_a; // extrinsic curvature

    arr_t & DIFFalpha_a = bssn->DIFFalpha_a;
    
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
          K_a(i, j, k) = K_src;
          DIFFchi_a(i,j,k) = DIFFchi[0][INDEX(i, j, k)];
          phi_a(i, j, k) = phi[INDEX(i, j, k)];
          DIFFalpha_a(i, j, k) = DIFFalpha_0;
        }
      }
    }

    for(int k = inner_lower[2]; k <= inner_upper[2]; k++)
    {
      for(int j = inner_lower[1]; j <= inner_upper[1]; j++)
      {
        for(int i = inner_lower[0]; i <= inner_upper[0]; i++)
        {
          psi1_a(i, j, k) = derivative(i, j, k, 1, phi_a, dx);
          psi2_a(i, j, k) = derivative(i, j, k, 2, phi_a, dx);
          psi3_a(i, j, k) = derivative(i, j, k, 3, phi_a, dx);
        }
      }
    }   
  }

  bssn->K0 = K_src;
  
  return flag;
}
   
}
