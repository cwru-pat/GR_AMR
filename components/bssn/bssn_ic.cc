#include "bssn_ic.h"
#include "../../cosmo_includes.h"

using namespace SAMRAI;

namespace cosmo
{

void bssn_ic_awa_stability(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln, real_t A)
{
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  std::random_device rd;
  std::mt19937 gen(7.0 /*rd()*/);
  
  boost::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
       hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;


     



  int rho = 0.02 / grid_geometry.getDx()[0];
  
  std::uniform_real_distribution<real_t> dist(
    -A/rho/rho, A/rho/rho
  );

  idx_t DIFFchi_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFchi"), variable_db->getContext("ACTIVE"));
  idx_t A11_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A11"), variable_db->getContext("ACTIVE"));
  idx_t A12_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A12"), variable_db->getContext("ACTIVE"));
  idx_t A13_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A13"), variable_db->getContext("ACTIVE"));
  idx_t A22_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A22"), variable_db->getContext("ACTIVE"));
  idx_t A23_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A23"), variable_db->getContext("ACTIVE"));
  idx_t A33_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A33"), variable_db->getContext("ACTIVE"));
  idx_t DIFFgamma11_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma11"), variable_db->getContext("ACTIVE"));
  idx_t DIFFgamma12_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma12"), variable_db->getContext("ACTIVE"));
  idx_t DIFFgamma13_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma13"), variable_db->getContext("ACTIVE"));
  idx_t DIFFgamma22_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma22"), variable_db->getContext("ACTIVE"));
  idx_t DIFFgamma23_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma23"), variable_db->getContext("ACTIVE"));
  idx_t DIFFgamma33_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma33"), variable_db->getContext("ACTIVE"));
  idx_t DIFFK_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFK"), variable_db->getContext("ACTIVE"));
  idx_t Gamma1_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("Gamma1"), variable_db->getContext("ACTIVE"));
  idx_t Gamma2_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("Gamma2"), variable_db->getContext("ACTIVE"));
  idx_t Gamma3_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("Gamma3"), variable_db->getContext("ACTIVE"));


  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;

    boost::shared_ptr<pdat::CellData<real_t> > DIFFchi_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFchi_a_idx)));

    boost::shared_ptr<pdat::CellData<real_t> > DIFFgamma11_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma11_a_idx)));
    boost::shared_ptr<pdat::CellData<real_t> > DIFFgamma12_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma12_a_idx)));
    boost::shared_ptr<pdat::CellData<real_t> > DIFFgamma13_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma13_a_idx)));
    boost::shared_ptr<pdat::CellData<real_t> > DIFFgamma22_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma22_a_idx)));
    boost::shared_ptr<pdat::CellData<real_t> > DIFFgamma23_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma23_a_idx)));
    boost::shared_ptr<pdat::CellData<real_t> > DIFFgamma33_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma33_a_idx)));

    boost::shared_ptr<pdat::CellData<real_t> > A11_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A11_a_idx)));
    boost::shared_ptr<pdat::CellData<real_t> > A12_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A12_a_idx)));
    boost::shared_ptr<pdat::CellData<real_t> > A13_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A13_a_idx)));
    boost::shared_ptr<pdat::CellData<real_t> > A22_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A22_a_idx)));
    boost::shared_ptr<pdat::CellData<real_t> > A23_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A23_a_idx)));
    boost::shared_ptr<pdat::CellData<real_t> > A33_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A33_a_idx)));

    boost::shared_ptr<pdat::CellData<real_t> > DIFFK_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFK_a_idx)));
    boost::shared_ptr<pdat::CellData<real_t> > Gamma1_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(Gamma1_a_idx)));
    boost::shared_ptr<pdat::CellData<real_t> > Gamma2_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(Gamma2_a_idx)));
    boost::shared_ptr<pdat::CellData<real_t> > Gamma3_a_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(Gamma3_a_idx)));

    arr_t DIFFchi_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFchi_a_pdata->getArrayData());     
    arr_t DIFFgamma11_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma11_a_pdata->getArrayData());     
    arr_t DIFFgamma12_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma12_a_pdata->getArrayData());     
    arr_t DIFFgamma13_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma13_a_pdata->getArrayData());     
    arr_t DIFFgamma22_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma22_a_pdata->getArrayData());     
    arr_t DIFFgamma23_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma23_a_pdata->getArrayData());     
    arr_t DIFFgamma33_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma33_a_pdata->getArrayData());     

    arr_t A11_a = pdat::ArrayDataAccess::access<DIM, double>(
      A11_a_pdata->getArrayData());     
    arr_t A12_a = pdat::ArrayDataAccess::access<DIM, double>(
      A12_a_pdata->getArrayData());     
    arr_t A13_a = pdat::ArrayDataAccess::access<DIM, double>(
      A13_a_pdata->getArrayData());     
    arr_t A22_a = pdat::ArrayDataAccess::access<DIM, double>(
      A22_a_pdata->getArrayData());     
    arr_t A23_a = pdat::ArrayDataAccess::access<DIM, double>(
      A23_a_pdata->getArrayData());     
    arr_t A33_a = pdat::ArrayDataAccess::access<DIM, double>(
      A33_a_pdata->getArrayData());     
    arr_t DIFFK_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFK_a_pdata->getArrayData());     
    arr_t Gamma1_a = pdat::ArrayDataAccess::access<DIM, double>(
      Gamma1_a_pdata->getArrayData());     
    arr_t Gamma2_a = pdat::ArrayDataAccess::access<DIM, double>(
      Gamma2_a_pdata->getArrayData());     
    arr_t Gamma3_a = pdat::ArrayDataAccess::access<DIM, double>(
      Gamma3_a_pdata->getArrayData());     

  
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

    
    const hier::Box& box = DIFFchi_a_pdata->getGhostBox();

    //const double *dx = &grid_geometry.getDx()[0];
      

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          DIFFgamma11_a(i,j,k) = dist(gen);
          DIFFgamma12_a(i,j,k) = dist(gen);
          DIFFgamma13_a(i,j,k) = dist(gen);
          DIFFgamma22_a(i,j,k) = dist(gen);
          DIFFgamma23_a(i,j,k) = dist(gen);
          DIFFgamma33_a(i,j,k) = dist(gen);
          DIFFchi_a(i,j,k) = dist(gen);
          A11_a(i,j,k)     = dist(gen);
          A12_a(i,j,k)     = dist(gen);
          A13_a(i,j,k)     = dist(gen);
          A22_a(i,j,k)     = dist(gen);
          A23_a(i,j,k)     = dist(gen);
          A33_a(i,j,k)     = dist(gen);
          DIFFK_a(i,j,k)   = dist(gen);
          Gamma1_a(i,j,k)  = dist(gen);
          Gamma2_a(i,j,k)  = dist(gen);
          Gamma3_a(i,j,k)  = dist(gen);
        }
      }
    }
                
  }
  
}

void bssn_ic_static_blackhole(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln)
{
# if ! USE_BSSN_SHIFT
  std::cerr << "Waning! USE_BSSN_SHIFT is suggested to be enabled in blackhole test " << std::endl;
# endif


  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

  idx_t DIFFchi_p_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFchi"), variable_db->getContext("PREVIOUS"));

  idx_t DIFFchi_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFchi"), variable_db->getContext("ACTIVE"));

  boost::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
       hierarchy->getGridGeometry()));
   TBOX_ASSERT(grid_geometry_);
   geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;


    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

    
    boost::shared_ptr<pdat::CellData<real_t> > chi_p_pdata(
      BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFchi_p_idx)));

    const hier::Box& box = chi_p_pdata->getGhostBox();
    
    boost::shared_ptr<pdat::CellData<real_t> > chi_a_pdata(
       BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFchi_a_idx)));

    // boost::shared_ptr<pdat::CellData<real_t> > chi_f_pdata(
    //   BOOST_CAST<pdat::CellData<real_t>, hier::PatchData>(
    //     patch->getPatchData(DIFFchi_f_idx)));


    arr_t chi_p =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        chi_p_pdata->getArrayData());

    arr_t chi_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        chi_a_pdata->getArrayData());

    const double * domain_lower = &grid_geometry.getXLower()[0];
    const double * domain_upper = &grid_geometry.getXUpper()[0];

    //const double *dx = &grid_geometry.getDx()[0];
    const double *dx = &patch_geom->getDx()[0];
      
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
          real_t x = (dx[0] * ((real_t)i + 0.5)) - L[0] / 2.0 ;
          real_t y = (dx[1] * ((real_t)j + 0.5)) - L[1] / 2.0 ;
          real_t z = (dx[2] * ((real_t)k + 0.5)) - L[2] / 2.0 ;

          real_t norm = sqrt(x*x + y*y + z*z);

          chi_p(i,j,k) = chi_a(i,j,k)
            //= 1/pw2((1.0 + 1.0/(2.0*norm))) - 1.0;
            = pw2(2.0*norm / (1+ 2.0 * norm)) - 1.0;
        }
      }
    }
                
  }
  
}
  
} // namespace cosmo
