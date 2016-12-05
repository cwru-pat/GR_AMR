#include "bssn_ic.h"
#include "../../cosmo_types.h"
#include "../../cosmo_globals.h"

namespace cosmo
{

/**
 * @brief      Set flat vacuum spacetime initial conditions,
 * plus small "noise" fluctuations.
 *
 * @param      bssn  BSSN instance
 */
void bssn_ic_awa_stability(BSSN * bssn, real_t A)
{
  idx_t i, j, k;

  std::random_device rd;
  std::mt19937 gen(7.0 /*rd()*/);
  std::uniform_real_distribution<real_t> dist(
      -A*50/N*50/N, A*50/N*50/N
    );

  arr_t & DIFFgamma11_p = *bssn->fields["DIFFgamma11_p"];
  arr_t & DIFFgamma12_p = *bssn->fields["DIFFgamma12_p"];
  arr_t & DIFFgamma13_p = *bssn->fields["DIFFgamma13_p"];
  arr_t & DIFFgamma22_p = *bssn->fields["DIFFgamma22_p"];
  arr_t & DIFFgamma23_p = *bssn->fields["DIFFgamma23_p"];
  arr_t & DIFFgamma33_p = *bssn->fields["DIFFgamma33_p"];
  arr_t & DIFFphi_p = *bssn->fields["DIFFphi_p"];
  arr_t & A11_p = *bssn->fields["A11_p"];
  arr_t & A12_p = *bssn->fields["A12_p"];
  arr_t & A13_p = *bssn->fields["A13_p"];
  arr_t & A22_p = *bssn->fields["A22_p"];
  arr_t & A23_p = *bssn->fields["A23_p"];
  arr_t & A33_p = *bssn->fields["A33_p"];
  arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];
  arr_t & Gamma1_p = *bssn->fields["Gamma1_p"];
  arr_t & Gamma2_p = *bssn->fields["Gamma2_p"];
  arr_t & Gamma3_p = *bssn->fields["Gamma3_p"];

  LOOP3(i, j, k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    DIFFgamma11_p[idx] = dist(gen);
    DIFFgamma12_p[idx] = dist(gen);
    DIFFgamma13_p[idx] = dist(gen);
    DIFFgamma22_p[idx] = dist(gen);
    DIFFgamma23_p[idx] = dist(gen);
    DIFFgamma33_p[idx] = dist(gen);
    DIFFphi_p[idx] = dist(gen);
    A11_p[idx]     = dist(gen);
    A12_p[idx]     = dist(gen);
    A13_p[idx]     = dist(gen);
    A22_p[idx]     = dist(gen);
    A23_p[idx]     = dist(gen);
    A33_p[idx]     = dist(gen);
    DIFFK_p[idx]   = dist(gen);
    Gamma1_p[idx]  = dist(gen);
    Gamma2_p[idx]  = dist(gen);
    Gamma3_p[idx]  = dist(gen);
  }
}

/**
 * @brief Initialize with a "linear" wave propagating in the x-direction.
 * @details AwA test.
 */
void bssn_ic_awa_linear_wave(BSSN * bssn)
{
  bssn_ic_awa_linear_wave(bssn, 1.0e-8, 1);
}

/**
 * @brief Initialize with a "linear" wave propagating in the 'dir'-direction.
 * @details AwA test.
 */
void bssn_ic_awa_linear_wave(BSSN * bssn, real_t A, int dir)
{
  idx_t i, j, k;

  arr_t & DIFFgamma11_p = *bssn->fields["DIFFgamma11_p"];
  arr_t & DIFFgamma22_p = *bssn->fields["DIFFgamma22_p"];
  arr_t & DIFFgamma33_p = *bssn->fields["DIFFgamma33_p"];

  arr_t & A11_p = *bssn->fields["A11_p"];
  arr_t & A22_p = *bssn->fields["A22_p"];
  arr_t & A33_p = *bssn->fields["A33_p"];

  LOOP3(i,j,k)
  {
    // position of dir
    real_t w = 0;
    switch(dir)
    {
      case 1 :
        w = ((real_t) i)*dx;
        DIFFgamma22_p[NP_INDEX(i,j,k)] = A*sin( 2.0*PI*w );
        DIFFgamma33_p[NP_INDEX(i,j,k)] = -A*sin( 2.0*PI*w );
        A22_p[NP_INDEX(i,j,k)] = PI*A*cos( 2.0*PI*w );
        A33_p[NP_INDEX(i,j,k)] = -PI*A*cos( 2.0*PI*w );
        break;
      case 2 :
        w = ((real_t) j)*dx;
        DIFFgamma11_p[NP_INDEX(i,j,k)] = A*sin( 2.0*PI*w );
        DIFFgamma33_p[NP_INDEX(i,j,k)] = -A*sin( 2.0*PI*w );
        A11_p[NP_INDEX(i,j,k)] = PI*A*cos( 2.0*PI*w );
        A33_p[NP_INDEX(i,j,k)] = -PI*A*cos( 2.0*PI*w );
        break;
      case 3 :
        w = ((real_t) k)*dx;
        DIFFgamma11_p[NP_INDEX(i,j,k)] = A*sin( 2.0*PI*w );
        DIFFgamma22_p[NP_INDEX(i,j,k)] = -A*sin( 2.0*PI*w );
        A11_p[NP_INDEX(i,j,k)] = PI*A*cos( 2.0*PI*w );
        A22_p[NP_INDEX(i,j,k)] = -PI*A*cos( 2.0*PI*w );
        break;
    }
  }
}

/**
 * @brief Initialize with a "linear" wave propagating in the x-direction.
 * @details AwA test.
 */
void bssn_ic_awa_diagonal_linear_wave(BSSN * bssn)
{
  // TODO
}

/**
 * @brief Initialize with a "linear" wave propagating in the x-direction,
 * on top of a deSitter spacetime.
 * @details Solution should behave, at linear order, as a damped wave equation.
 */
void bssn_ic_awa_linear_wave_desitter(BSSN * bssn)
{
  bssn_ic_awa_linear_wave(bssn);

  idx_t i, j, k;

  arr_t & r_a = *bssn->fields["r_a"];
  arr_t & K_p = *bssn->fields["K_p"];

  LOOP3(i,j,k)
  {
    // FRW Background parameters
    real_t rho = 1.0;

    r_a[NP_INDEX(i,j,k)] = rho; // constant density
    K_p[NP_INDEX(i,j,k)] = -sqrt(24.0*PI*rho);
  }
}

/**
 * @brief      AwA Gauge Wave Test
 */
void bssn_ic_awa_gauge_wave(BSSN * bssn)
{
  bssn_ic_awa_gauge_wave(bssn, 1);
}

/**
 * @brief      AwA Gauge Wave Test, setting wave propagation direction
 */
void bssn_ic_awa_gauge_wave(BSSN * bssn, int dir)
{
  idx_t i, j, k;

  real_t A = 0.5;

  arr_t & DIFFphi_p = *bssn->fields["DIFFphi_p"];
  arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];
  arr_t & DIFFalpha_p = *bssn->fields["DIFFalpha_p"];

  arr_t & DIFFgamma11_p = *bssn->fields["DIFFgamma11_p"];
  arr_t & DIFFgamma22_p = *bssn->fields["DIFFgamma22_p"];
  arr_t & DIFFgamma33_p = *bssn->fields["DIFFgamma33_p"];
  
  arr_t & A11_p = *bssn->fields["A11_p"];
  arr_t & A22_p = *bssn->fields["A22_p"];
  arr_t & A33_p = *bssn->fields["A33_p"];

  arr_t & Gamma1_p = *bssn->fields["Gamma1_p"];
  arr_t & Gamma2_p = *bssn->fields["Gamma2_p"];
  arr_t & Gamma3_p = *bssn->fields["Gamma3_p"];

  LOOP3(i,j,k)
  {
    // position of dir
    real_t w = 0;
    switch(dir)
    {
      case 1 :
        w = ((real_t) i)*dx;
        break;
      case 2 :
        w = ((real_t) j)*dx;
        break;
      case 3 :
        w = ((real_t) k)*dx;
        break;
    }

    // H(t = 0)
    real_t H = A*sin( 2.0*PI*w );
    real_t dtH = -2.0*PI*A*cos( 2.0*PI*w );
    real_t Kxx = dtH/(2.0*sqrt(1.0-H));
    real_t K = Kxx / (1.0 - H);

    DIFFphi_p[NP_INDEX(i,j,k)] = log1p(-H)/12.0;
    DIFFalpha_p[NP_INDEX(i,j,k)] = sqrt(1.0-H) - 1.0;
    DIFFK_p[NP_INDEX(i,j,k)] = K;

    DIFFgamma11_p[NP_INDEX(i,j,k)] = pow(1.0-H, -1.0/3.0) - 1.0;
    DIFFgamma22_p[NP_INDEX(i,j,k)] = pow(1.0-H, -1.0/3.0) - 1.0;
    DIFFgamma33_p[NP_INDEX(i,j,k)] = pow(1.0-H, -1.0/3.0) - 1.0;

    A11_p[NP_INDEX(i,j,k)] = - K/3.0*pow(1.0-H, -1.0/3.0);
    A22_p[NP_INDEX(i,j,k)] = - K/3.0*pow(1.0-H, -1.0/3.0);
    A33_p[NP_INDEX(i,j,k)] = - K/3.0*pow(1.0-H, -1.0/3.0);

    switch(dir)
    {
      case 1 :
        DIFFgamma11_p[NP_INDEX(i,j,k)] = pow(1.0-H, 2.0/3.0) - 1.0;
        A11_p[NP_INDEX(i,j,k)] = Kxx*pow(1.0-H, -1.0/3.0) - K/3.0*pow(1.0-H, 2.0/3.0);
        Gamma1_p[NP_INDEX(i,j,k)] = 2.0/3.0*dtH*pow(1.0-H, -5.0/3.0);
        break;
      case 2 :
        DIFFgamma22_p[NP_INDEX(i,j,k)] = pow(1.0-H, 2.0/3.0) - 1.0;
        A22_p[NP_INDEX(i,j,k)] = Kxx*pow(1.0-H, -1.0/3.0) - K/3.0*pow(1.0-H, 2.0/3.0);
        Gamma2_p[NP_INDEX(i,j,k)] = 2.0/3.0*dtH*pow(1.0-H, -5.0/3.0);
        break;
      case 3 :
        DIFFgamma33_p[NP_INDEX(i,j,k)] = pow(1.0-H, 2.0/3.0) - 1.0;
        A33_p[NP_INDEX(i,j,k)] = Kxx*pow(1.0-H, -1.0/3.0) - K/3.0*pow(1.0-H, 2.0/3.0);
        Gamma3_p[NP_INDEX(i,j,k)] = 2.0/3.0*dtH*pow(1.0-H, -5.0/3.0);
        break;
      default:
        throw -1;
        break;
    }

  }
}

void bssn_ic_awa_shifted_gauge_wave(BSSN * bssn)
{
  bssn_ic_awa_shifted_gauge_wave(bssn, 1);
}

/**
 * @brief      AwA Shifted Gauge Wave Test
 */
void bssn_ic_awa_shifted_gauge_wave(BSSN * bssn, int dir)
{
# if ! USE_BSSN_SHIFT
  std::cerr << "USE_BSSN_SHIFT must be enabled for shifted gauge wave! (cmake using -DCOSMO_USE_BSSN_SHIFT=1)" << std::endl;
  throw -1;
# endif

  idx_t i, j, k;

  real_t A = 0.5;

  arr_t & DIFFphi_p = *bssn->fields["DIFFphi_p"];
  arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];
  arr_t & DIFFalpha_p = *bssn->fields["DIFFalpha_p"];

  arr_t & DIFFgamma11_p = *bssn->fields["DIFFgamma11_p"];
  arr_t & DIFFgamma22_p = *bssn->fields["DIFFgamma22_p"];
  arr_t & DIFFgamma33_p = *bssn->fields["DIFFgamma33_p"];
  
  arr_t & A11_p = *bssn->fields["A11_p"];
  arr_t & A22_p = *bssn->fields["A22_p"];
  arr_t & A33_p = *bssn->fields["A33_p"];

  arr_t & Gamma1_p = *bssn->fields["Gamma1_p"];
  arr_t & Gamma2_p = *bssn->fields["Gamma2_p"];
  arr_t & Gamma3_p = *bssn->fields["Gamma3_p"];

  arr_t & beta1_p = *bssn->fields["beta1_p"];
  arr_t & beta2_p = *bssn->fields["beta2_p"];
  arr_t & beta3_p = *bssn->fields["beta3_p"];

  LOOP3(i,j,k)
  {
    // position of dir
    real_t w = 0;
    switch(dir)
    {
      case 1 :
        w = ((real_t) i)*dx;
        break;
      case 2 :
        w = ((real_t) j)*dx;
        break;
      case 3 :
        w = ((real_t) k)*dx;
        break;
    }

    // H(t = 0)
    real_t H = A*sin( 2.0*PI*w );
    real_t dtH = -2.0*PI*A*cos( 2.0*PI*w );
    real_t Kxx = dtH/(2.0*sqrt(1.0+H));
    real_t K = Kxx / (1.0 + H);

    DIFFphi_p[NP_INDEX(i,j,k)] = log1p(H)/12.0;
    DIFFalpha_p[NP_INDEX(i,j,k)] = sqrt(1.0/(1.0+H)) - 1.0;
    DIFFK_p[NP_INDEX(i,j,k)] = K;

    DIFFgamma11_p[NP_INDEX(i,j,k)] = pow(1.0+H, -1.0/3.0) - 1.0;
    DIFFgamma22_p[NP_INDEX(i,j,k)] = pow(1.0+H, -1.0/3.0) - 1.0;
    DIFFgamma33_p[NP_INDEX(i,j,k)] = pow(1.0+H, -1.0/3.0) - 1.0;

    A11_p[NP_INDEX(i,j,k)] = - K/3.0*pow(1.0+H, -1.0/3.0);
    A22_p[NP_INDEX(i,j,k)] = - K/3.0*pow(1.0+H, -1.0/3.0);
    A33_p[NP_INDEX(i,j,k)] = - K/3.0*pow(1.0+H, -1.0/3.0);

    switch(dir)
    {
      case 1 :
        Gamma1_p[NP_INDEX(i,j,k)] = -2.0/3.0*dtH*pow(1.0+H, -5.0/3.0);
        beta1_p[NP_INDEX(i,j,k)] = -H/(1.0+H);
        DIFFgamma11_p[NP_INDEX(i,j,k)] = pow(1.0+H, 2.0/3.0) - 1.0;
        A11_p[NP_INDEX(i,j,k)] = Kxx*pow(1.0+H, -1.0/3.0) - K/3.0*pow(1.0+H, 2.0/3.0);
        break;
      case 2 :
        Gamma2_p[NP_INDEX(i,j,k)] = -2.0/3.0*dtH*pow(1.0+H, -5.0/3.0);
        beta2_p[NP_INDEX(i,j,k)] = -H/(1.0+H);
        DIFFgamma22_p[NP_INDEX(i,j,k)] = pow(1.0+H, 2.0/3.0) - 1.0;
        A22_p[NP_INDEX(i,j,k)] = Kxx*pow(1.0+H, -1.0/3.0) - K/3.0*pow(1.0+H, 2.0/3.0);
        break;
      case 3 :
        Gamma3_p[NP_INDEX(i,j,k)] = -2.0/3.0*dtH*pow(1.0+H, -5.0/3.0);
        beta3_p[NP_INDEX(i,j,k)] = -H/(1.0+H);
        DIFFgamma33_p[NP_INDEX(i,j,k)] = pow(1.0+H, 2.0/3.0) - 1.0;
        A33_p[NP_INDEX(i,j,k)] = Kxx*pow(1.0+H, -1.0/3.0) - K/3.0*pow(1.0+H, 2.0/3.0);
        break;
      default:
        throw -1;
        break;
    }
  }
}

void bssn_ic_static_blackhole(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
# if ! USE_BSSN_SHIFT
  std::cerr << "Waning! USE_BSSN_SHIFT is suggested to be enabled in blackhole test " << std::endl;
# endif

  idx_t i, j, k;

  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

  idx_t chi_p_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFchi"), variable_db->getContext("PREVIOUS"));


  boost::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(0));

  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
       hierarchy->getGridGeometry()));
   TBOX_ASSERT(grid_geometry_);
   geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;

    const hier::Box& box = patch->getBox();

    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));
    
    boost::shared_ptr<pdat::CellData<real_t> > chi_p_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFchi_p_idx)));

    boost::shared_ptr<pdat::CellData<real_t> > chi_a_pdata(
       BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFchi_a_idx)));

    // boost::shared_ptr<pdat::CellData<real_t> > chi_f_pdata(
    //   BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
    //     patch->getPatchData(DIFFchi_f_idx)));


    arr_t chi_p =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        chi_p_pdata->getArrayData());
    
    const double * domain_lower = &grid_geometry.getXLower()[0];
    const double * domain_upper = &grid_geometry.getXUpper()[0];

    const double *dx = &grid_geometry.getDx()[0];

      
    double L[3];


    for(int i = 0 ; i < 3; i++)
      L[i] = domain_upper[i] - domain_lower[i];

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          real_t x = (dx[0] * (real_t)(((i+(N-INNER_N)/2)%N)+0.5)) - L[0] / 2.0 ;
          real_t y = (dx[1] * (real_t)(((j+(N-INNER_N)/2)%N)+0.5)) - L[1] / 2.0 ;
          real_t z = (dx[2] * (real_t)(((k+(N-INNER_N)/2)%N)+0.5)) - L[2] / 2.0 ;

          real_t norm = sqrt(x*x + y*y + z*z);

          chi_p(i,j,k) = log(1.0 + 1.0/(2.0*norm));

        }
      }
    }
    
            
  }
  




  
  LOOP3(i,j,k)
  {
    real_t x = (dx * (real_t)(((i+(N-INNER_N)/2)%N)+0.5)) - H_LEN_FRAC / 2.0 ;
    real_t y = (dx * (real_t)(((j+(N-INNER_N)/2)%N)+0.5)) - H_LEN_FRAC / 2.0 ;
    real_t z = (dx * (real_t)(((k+(N-INNER_N)/2)%N)+0.5)) - H_LEN_FRAC / 2.0 ;

    real_t norm = sqrt(x*x + y*y + z*z);

    DIFFphi_p[NP_INDEX(i,j,k)] = log(1.0 + 1.0/(2.0*norm));
  }
}
  
} // namespace cosmo
