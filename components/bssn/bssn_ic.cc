#include "bssn_ic.h"
#include "../../cosmo_includes.h"
#include "SAMRAI/tbox/HDFDatabase.h"
#include "../elliptic_solver/full_multigrid.h"
#include "../elliptic_solver/multigrid_bd_handler.h"

using namespace SAMRAI;

namespace cosmo
{

void bssn_ic_awa_stability(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln, real_t A)
{
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  std::random_device rd;
  std::mt19937 gen(7.0 /*rd()*/);
  
  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
       hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;


     



  double rho = 0.02 / grid_geometry.getDx()[0];
  
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
    const std::shared_ptr<hier::Patch> & patch = *pit;

    std::shared_ptr<pdat::CellData<real_t> > DIFFchi_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFchi_a_idx)));

    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma12_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma12_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma13_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma13_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma22_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma23_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma23_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma33_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma33_a_idx)));

    std::shared_ptr<pdat::CellData<real_t> > A11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A12_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A12_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A13_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A13_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A22_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A23_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A23_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A33_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A33_a_idx)));

    std::shared_ptr<pdat::CellData<real_t> > DIFFK_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFK_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > Gamma1_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(Gamma1_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > Gamma2_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(Gamma2_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > Gamma3_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
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

  
    const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

    
    const hier::Box& box = DIFFchi_a_pdata->getBox();

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
  
void bssn_ic_awa_linear_wave(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln, real_t A, idx_t dir)
{
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  idx_t A11_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A11"), variable_db->getContext("ACTIVE"));
  idx_t A22_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A22"), variable_db->getContext("ACTIVE"));
  idx_t A33_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A33"), variable_db->getContext("ACTIVE"));

  idx_t DIFFgamma11_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma11"), variable_db->getContext("ACTIVE"));
  idx_t DIFFgamma22_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma22"), variable_db->getContext("ACTIVE"));
  idx_t DIFFgamma33_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma33"), variable_db->getContext("ACTIVE"));

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    std::shared_ptr<pdat::CellData<real_t> > A11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A22_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A33_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A33_a_idx)));

    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma22_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma33_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma33_a_idx)));

    arr_t DIFFgamma11_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma11_a_pdata->getArrayData());     
    arr_t DIFFgamma22_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma22_a_pdata->getArrayData());     
    arr_t DIFFgamma33_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma33_a_pdata->getArrayData());     
    arr_t A11_a = pdat::ArrayDataAccess::access<DIM, double>(
      A11_a_pdata->getArrayData());     
    arr_t A22_a = pdat::ArrayDataAccess::access<DIM, double>(
      A22_a_pdata->getArrayData());     
    arr_t A33_a = pdat::ArrayDataAccess::access<DIM, double>(
      A33_a_pdata->getArrayData());     

    const hier::Box& box = patch->getBox();

    const double *dx = &grid_geometry.getDx()[0];
    
    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
              real_t w = 0;
              switch(dir)
              {
              case 1 :
                w = ((real_t) i + 0.5)*dx[0] ;
                DIFFgamma22_a(i,j,k) = A*sin( 2.0*PI*w );
                DIFFgamma33_a(i,j,k) = -A*sin( 2.0*PI*w );
                A22_a(i,j,k) = PI*A*cos( 2.0*PI*w );
                A33_a(i,j,k) = -PI*A*cos( 2.0*PI*w );
                break;
              case 2 :
                w = ((real_t) j + 0.5)*dx[1];
                DIFFgamma11_a(i,j,k) = A*sin( 2.0*PI*w );
                DIFFgamma33_a(i,j,k) = -A*sin( 2.0*PI*w );
                A11_a(i,j,k) = PI*A*cos( 2.0*PI*w );
                A33_a(i,j,k) = -PI*A*cos( 2.0*PI*w );
                break;
              case 3 :
                w = ((real_t) k + 0.5)*dx[2];
                DIFFgamma11_a(i,j,k) = A*sin( 2.0*PI*w );
                DIFFgamma22_a(i,j,k) = -A*sin( 2.0*PI*w );
                A11_a(i,j,k) = PI*A*cos( 2.0*PI*w );
                A22_a(i,j,k) = -PI*A*cos( 2.0*PI*w );
                break;
              }

        }
      }
    }


  }
}

/**
 * @brief      AwA Gauge Wave Test, setting wave propagation direction
 */
void bssn_ic_awa_gauge_wave(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln, idx_t dir)
{
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  real_t A = 0.5;
  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  idx_t DIFFchi_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFchi"), variable_db->getContext("ACTIVE"));
  idx_t DIFFK_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFK"), variable_db->getContext("ACTIVE"));
  idx_t DIFFalpha_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFalpha"), variable_db->getContext("ACTIVE"));

  idx_t A11_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A11"), variable_db->getContext("ACTIVE"));
  idx_t A22_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A22"), variable_db->getContext("ACTIVE"));
  idx_t A33_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A33"), variable_db->getContext("ACTIVE"));
  idx_t DIFFgamma11_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma11"), variable_db->getContext("ACTIVE"));
  idx_t DIFFgamma22_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma22"), variable_db->getContext("ACTIVE"));
  idx_t DIFFgamma33_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma33"), variable_db->getContext("ACTIVE"));
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
    const std::shared_ptr<hier::Patch> & patch = *pit;

    std::shared_ptr<pdat::CellData<real_t> > DIFFchi_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFchi_a_idx)));

    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma22_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma33_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma33_a_idx)));
    
    std::shared_ptr<pdat::CellData<real_t> > A11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A22_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A33_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A33_a_idx)));

    std::shared_ptr<pdat::CellData<real_t> > DIFFK_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFK_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > Gamma1_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(Gamma1_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > Gamma2_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(Gamma2_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > Gamma3_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(Gamma3_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFalpha_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFalpha_a_idx)));

    arr_t DIFFchi_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFchi_a_pdata->getArrayData());     
    arr_t DIFFgamma11_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma11_a_pdata->getArrayData());
    arr_t DIFFgamma22_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma22_a_pdata->getArrayData());     
        
    arr_t DIFFgamma33_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma33_a_pdata->getArrayData());     

    arr_t A11_a = pdat::ArrayDataAccess::access<DIM, double>(
      A11_a_pdata->getArrayData());     
    arr_t A22_a = pdat::ArrayDataAccess::access<DIM, double>(
      A22_a_pdata->getArrayData());     
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
    arr_t DIFFalpha_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFalpha_a_pdata->getArrayData());     

  
    const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry()));

    
    const hier::Box& box = DIFFchi_a_pdata->getBox();
    const double *dx = &grid_geometry.getDx()[0];
      

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          // position of dir
          real_t w = 0;
          switch(dir)
          {
          case 1 :
            w = ((real_t) i + 0.5)*dx[0];
            break;
          case 2 :
            w = ((real_t) j + 0.5)*dx[1];
            break;
          case 3 :
            w = ((real_t) k + 0.5)*dx[2];
            break;
          }

          // H(t = 0)
          real_t H = A*sin( 2.0*PI*w );
          real_t dtH = -2.0*PI*A*cos( 2.0*PI*w );
          real_t Kxx = dtH/(2.0*sqrt(1.0-H));
          real_t K = Kxx / (1.0 - H);

          DIFFchi_a(i,j,k) = pow(-H+1.0, -1.0/6.0) - 1.0;
          DIFFalpha_a(i,j,k) = sqrt(1.0-H) - 1.0;
          DIFFK_a(i,j,k) = K;

          DIFFgamma11_a(i,j,k) = pow(1.0-H, -1.0/3.0) - 1.0;
          DIFFgamma22_a(i,j,k) = pow(1.0-H, -1.0/3.0) - 1.0;
          DIFFgamma33_a(i,j,k) = pow(1.0-H, -1.0/3.0) - 1.0;

          A11_a(i,j,k) = - K/3.0*pow(1.0-H, -1.0/3.0);
          A22_a(i,j,k) = - K/3.0*pow(1.0-H, -1.0/3.0);
          A33_a(i,j,k) = - K/3.0*pow(1.0-H, -1.0/3.0);

          switch(dir)
          {
          case 1 :
            DIFFgamma11_a(i,j,k) = pow(1.0-H, 2.0/3.0) - 1.0;
            A11_a(i,j,k) = Kxx*pow(1.0-H, -1.0/3.0) - K/3.0*pow(1.0-H, 2.0/3.0);
            Gamma1_a(i,j,k) = 2.0/3.0*dtH*pow(1.0-H, -5.0/3.0);
            break;
          case 2 :
            DIFFgamma22_a(i,j,k) = pow(1.0-H, 2.0/3.0) - 1.0;
            A22_a(i,j,k) = Kxx*pow(1.0-H, -1.0/3.0) - K/3.0*pow(1.0-H, 2.0/3.0);
            Gamma2_a(i,j,k) = 2.0/3.0*dtH*pow(1.0-H, -5.0/3.0);
            break;
          case 3 :
            DIFFgamma33_a(i,j,k) = pow(1.0-H, 2.0/3.0) - 1.0;
            A33_a(i,j,k) = Kxx*pow(1.0-H, -1.0/3.0) - K/3.0*pow(1.0-H, 2.0/3.0);
            Gamma3_a(i,j,k) = 2.0/3.0*dtH*pow(1.0-H, -5.0/3.0);
            break;
          default:
            throw -1;
            break;
      
          }
        }
      }
    }
  }
}

/**
 * @brief      AwA Gauge Wave Test, setting wave propagation direction
 */
void bssn_ic_awa_shifted_gauge_wave(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln, idx_t dir)
{
# if ! USE_BSSN_SHIFT
  TBOX_ERROR("USE_BSSN_SHIFT must be enabled for shifted gauge wave! (cmake using -DCOSMO_USE_BSSN_SHIFT=1)\n");
# endif
  
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  real_t A = 0.5;
  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  idx_t DIFFchi_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFchi"), variable_db->getContext("ACTIVE"));
  idx_t DIFFK_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFK"), variable_db->getContext("ACTIVE"));
  idx_t DIFFalpha_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFalpha"), variable_db->getContext("ACTIVE"));

  idx_t A11_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A11"), variable_db->getContext("ACTIVE"));
  idx_t A22_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A22"), variable_db->getContext("ACTIVE"));
  idx_t A33_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("A33"), variable_db->getContext("ACTIVE"));
  idx_t DIFFgamma11_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma11"), variable_db->getContext("ACTIVE"));
  idx_t DIFFgamma22_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma22"), variable_db->getContext("ACTIVE"));
  idx_t DIFFgamma33_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFgamma33"), variable_db->getContext("ACTIVE"));
  idx_t Gamma1_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("Gamma1"), variable_db->getContext("ACTIVE"));
  idx_t Gamma2_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("Gamma2"), variable_db->getContext("ACTIVE"));
  idx_t Gamma3_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("Gamma3"), variable_db->getContext("ACTIVE"));
  idx_t beta1_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("beta1"), variable_db->getContext("ACTIVE"));
  idx_t beta2_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("beta2"), variable_db->getContext("ACTIVE"));
  idx_t beta3_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("beta3"), variable_db->getContext("ACTIVE"));

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    std::shared_ptr<pdat::CellData<real_t> > DIFFchi_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFchi_a_idx)));

    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma22_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma33_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma33_a_idx)));
    
    std::shared_ptr<pdat::CellData<real_t> > A11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A22_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A33_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A33_a_idx)));

    std::shared_ptr<pdat::CellData<real_t> > DIFFK_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFK_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > Gamma1_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(Gamma1_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > Gamma2_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(Gamma2_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > Gamma3_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(Gamma3_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > beta1_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(beta1_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > beta2_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(beta2_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > beta3_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(beta3_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFalpha_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFalpha_a_idx)));

    arr_t DIFFchi_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFchi_a_pdata->getArrayData());     
    arr_t DIFFgamma11_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma11_a_pdata->getArrayData());
    arr_t DIFFgamma22_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma22_a_pdata->getArrayData());     
        
    arr_t DIFFgamma33_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFgamma33_a_pdata->getArrayData());     

    arr_t A11_a = pdat::ArrayDataAccess::access<DIM, double>(
      A11_a_pdata->getArrayData());     
    arr_t A22_a = pdat::ArrayDataAccess::access<DIM, double>(
      A22_a_pdata->getArrayData());     
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
    arr_t DIFFalpha_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFalpha_a_pdata->getArrayData());     
    arr_t beta1_a = pdat::ArrayDataAccess::access<DIM, double>(
      beta1_a_pdata->getArrayData());     
    arr_t beta2_a = pdat::ArrayDataAccess::access<DIM, double>(
      beta2_a_pdata->getArrayData());     
    arr_t beta3_a = pdat::ArrayDataAccess::access<DIM, double>(
      beta3_a_pdata->getArrayData());     

  
    const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry()));

    
    const hier::Box& box = DIFFchi_a_pdata->getBox();
    const double *dx = &grid_geometry.getDx()[0];
      

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          // position of dir
          real_t w = 0;
          switch(dir)
          {
          case 1 :
            w = ((real_t) i + 0.5)*dx[0];
            break;
          case 2 :
            w = ((real_t) j + 0.5)*dx[1];
            break;
          case 3 :
            w = ((real_t) k + 0.5)*dx[2];
            break;
          }

          // H(t = 0)
          real_t H = A*sin( 2.0*PI*w );
          real_t dtH = -2.0*PI*A*cos( 2.0*PI*w );
          real_t Kxx = dtH/(2.0*sqrt(1.0+H));
          real_t K = Kxx / (1.0 + H);

          DIFFchi_a(i,j,k) = pow(H+1.0, -1.0/6.0) - 1.0;
          DIFFalpha_a(i,j,k) = sqrt(1.0/(1.0+H)) - 1.0;
          DIFFK_a(i,j,k) = K;

          DIFFgamma11_a(i,j,k) = pow(1.0+H, -1.0/3.0) - 1.0;
          DIFFgamma22_a(i,j,k) = pow(1.0+H, -1.0/3.0) - 1.0;
          DIFFgamma33_a(i,j,k) = pow(1.0+H, -1.0/3.0) - 1.0;

          A11_a(i,j,k) = - K/3.0*pow(1.0+H, -1.0/3.0);
          A22_a(i,j,k) = - K/3.0*pow(1.0+H, -1.0/3.0);
          A33_a(i,j,k) = - K/3.0*pow(1.0+H, -1.0/3.0);

          switch(dir)
          {
          case 1 :
            DIFFgamma11_a(i,j,k) = pow(1.0+H, 2.0/3.0) - 1.0;
            A11_a(i,j,k) = Kxx*pow(1.0+H, -1.0/3.0) - K/3.0*pow(1.0+H, 2.0/3.0);
            beta1_a(i,j,k) = -H/(1.0+H);
            Gamma1_a(i,j,k) = -2.0/3.0*dtH*pow(1.0+H, -5.0/3.0);
            break;
          case 2 :
            DIFFgamma22_a(i,j,k) = pow(1.0+H, 2.0/3.0) - 1.0;
            A22_a(i,j,k) = Kxx*pow(1.0+H, -1.0/3.0) - K/3.0*pow(1.0+H, 2.0/3.0);
            beta2_a(i,j,k) = -H/(1.0+H);
            Gamma2_a(i,j,k) = -2.0/3.0*dtH*pow(1.0+H, -5.0/3.0);
            break;
          case 3 :
            DIFFgamma33_a(i,j,k) = pow(1.0+H, 2.0/3.0) - 1.0;
            A33_a(i,j,k) = Kxx*pow(1.0+H, -1.0/3.0) - K/3.0*pow(1.0+H, 2.0/3.0);
            beta3_a(i,j,k) = -H/(1.0+H);
            Gamma3_a(i,j,k) = -2.0/3.0*dtH*pow(1.0+H, -5.0/3.0);
            break;
          default:
            throw -1;
            break;
      
          }
        }
      }
    }
  }
}

void bssn_ic_static_BHL_CTT(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln,
  real_t M,
  real_t a,
  real_t K_c,
  real_t relaxation_tolerance,
  idx_t num_vcycles)
{
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  idx_t DIFFchi_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFchi"), variable_db->getContext("ACTIVE"));

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

  idx_t DIFFK_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFK"), variable_db->getContext("ACTIVE"));

  
  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
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


   for(int i = 0 ; i < 3; i++)
     L[i] = domain_upper[i] - domain_lower[i];

   double l = L[0]/2 - 4.0* M;
   double sigma = 3.5 * M;

   std::string boundary_type = "periodic";
   multigridBdHandler * bd_handler = new multigridBdHandler(boundary_type, L, 10);

    
   idx_t NX = round(L[0] / dx[0]);
   idx_t NY = round(L[1] / dx[1]); 
   idx_t NZ = round(L[2] / dx[2]);


    // initializing solution vector
   CosmoArray<idx_t, real_t> * X = new CosmoArray<idx_t, real_t> [4];

   bool flag = false;

    X[0].init(NX, NY, NZ);
    X[1].init(NX, NY, NZ);
    X[2].init(NX, NY, NZ);
    X[3].init(NX, NY, NZ);

    idx_t molecule_n[4] = {18, 5, 5, 5};

    atom atom_tmp = {0};
    
    FASMultigrid multigrid(
      X, 4, molecule_n, 4, 2, relaxation_tolerance, L, NX, NY, NZ, bd_handler);

    /***********start initializaing equations **********************/

    multigrid.eqns[0][0].init(1, 1);
    //adding laplacian term
    atom_tmp.type =  multigrid.atom_type::lap;
    atom_tmp.u_id = 0;
    multigrid.eqns[0][0].add_atom(atom_tmp);

    //adding \Psi^-7 term
    for(int i = 1; i <= 3; i++)
      for(int j = 1; j <= 3; j++)
      {
        atom_tmp.type = i+1;
        atom_tmp.u_id = j;
        if(i == j)
          multigrid.eqns[0][3*(i-1)+j].init(3, 0.5 - 1.0/6.0);
        else
          multigrid.eqns[0][3*(i-1)+j].init(3, 0.25);
        multigrid.eqns[0][3*(i-1)+j].add_atom(atom_tmp);

        multigrid.eqns[0][3*(i-1)+j].add_atom(atom_tmp);

        atom_tmp.type = 1;
        atom_tmp.u_id = 0;
        atom_tmp.value = -7;
        multigrid.eqns[0][3*(i-1)+j].add_atom(atom_tmp);
      }
    for(int i = 1; i <= 3; i++)
      for(int j = i+1; j <= 3; j++)
      {
        atom_tmp.type = i+1;
        atom_tmp.u_id = j;
        multigrid.eqns[0][(9 + i+j - 2)].init(3, 0.5);
        multigrid.eqns[0][(9 + i+j - 2)].add_atom(atom_tmp);
      
        atom_tmp.type = j+1;
        atom_tmp.u_id = i;
        multigrid.eqns[0][(9 + i+j - 2)].add_atom(atom_tmp);

        atom_tmp.type = 1;
        atom_tmp.u_id = 0;
        atom_tmp.value = -7;
        multigrid.eqns[0][(9 + i+j - 2)].add_atom(atom_tmp);
      }

    for(int i = 1; i <= 3; i++)
      for(int j = i+1; j <= 3; j++)
      {
        multigrid.eqns[0][(12 + i+j - 2)].init(3, -1.0/3.0);
        atom_tmp.type = i+1;
        atom_tmp.u_id = i;
        multigrid.eqns[0][(12+i+j-2)].add_atom(atom_tmp);

        atom_tmp.type = j+1;
        atom_tmp.u_id = j;
        multigrid.eqns[0][(12+i+j-2)].add_atom(atom_tmp);

        atom_tmp.type = 1;
        atom_tmp.u_id = 0;
        atom_tmp.value = -7;
        multigrid.eqns[0][(12 + i+j - 2)].add_atom(atom_tmp);
      }

  
    multigrid.eqns[0][16].init(1, 1.0);

    //    atom_tmp.type = multigrid.atom_type::const_f;
    //multigrid.eqns[0][16].add_atom(atom_tmp);
  
    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 5;
    multigrid.eqns[0][16].add_atom(atom_tmp);

    multigrid.eqns[0][17].init(1, 1.0);
    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 0;
    multigrid.eqns[0][17].add_atom(atom_tmp);
    
    // start adding momentum constraints
    multigrid.eqns[1][0].init(1, 1.0);
    multigrid.eqns[1][1].init(1, 1.0/3.0);
    multigrid.eqns[1][2].init(1, 1.0/3.0);
    multigrid.eqns[1][3].init(1, 1.0/3.0);
    multigrid.eqns[1][4].init(1, 1.0);

    multigrid.eqns[2][0].init(1, 1.0);
    multigrid.eqns[2][1].init(1, 1.0/3.0);
    multigrid.eqns[2][2].init(1, 1.0/3.0);
    multigrid.eqns[2][3].init(1, 1.0/3.0);
    multigrid.eqns[2][4].init(1, 1.0);

    multigrid.eqns[3][0].init(1, 1.0);
    multigrid.eqns[3][1].init(1, 1.0/3.0);
    multigrid.eqns[3][2].init(1, 1.0/3.0);
    multigrid.eqns[3][3].init(1, 1.0/3.0);
    multigrid.eqns[3][4].init(1, 1.0);

    //adding terms to eqn 1
    atom_tmp.type = multigrid.atom_type::lap;
    atom_tmp.u_id = 1;
    multigrid.eqns[1][0].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der11;
    atom_tmp.u_id = 1;
    multigrid.eqns[1][1].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der12;
    atom_tmp.u_id = 2;
    multigrid.eqns[1][2].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der13;
    atom_tmp.u_id = 3;
    multigrid.eqns[1][3].add_atom(atom_tmp);

    //    atom_tmp.type = multigrid.atom_type::const_f;
    //multigrid.eqns[1][4].add_atom(atom_tmp);
    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 6;
    multigrid.eqns[1][4].add_atom(atom_tmp);

  
    //adding terms to eqn 2
    atom_tmp.type = multigrid.atom_type::lap;
    atom_tmp.u_id = 2;
    multigrid.eqns[2][0].add_atom(atom_tmp);


    atom_tmp.type = multigrid.atom_type::der12;
    atom_tmp.u_id = 1;
    multigrid.eqns[2][1].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der22;
    atom_tmp.u_id = 2;
    multigrid.eqns[2][2].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der23;
    atom_tmp.u_id = 3;
    multigrid.eqns[2][3].add_atom(atom_tmp);

    //    atom_tmp.type = multigrid.atom_type::const_f;
    //multigrid.eqns[2][4].add_atom(atom_tmp);
    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 6;
    multigrid.eqns[2][4].add_atom(atom_tmp);
 
    //adding terms to eqn 3
    atom_tmp.type = multigrid.atom_type::lap;
    atom_tmp.u_id = 3;
    multigrid.eqns[3][0].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der13;
    atom_tmp.u_id = 1;
    multigrid.eqns[3][1].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der23;
    atom_tmp.u_id = 2;
    multigrid.eqns[3][2].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der33;
    atom_tmp.u_id = 3;
    multigrid.eqns[3][3].add_atom(atom_tmp);

    //    atom_tmp.type = multigrid.atom_type::const_f;
    //multigrid.eqns[3][4].add_atom(atom_tmp);
    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 6;
    multigrid.eqns[3][4].add_atom(atom_tmp);

  
  /*********ending initializing equations************/
    
    for(int i=0; i<NX; ++i) 
    for(int j=0; j<NY; ++j) 
    for(int k=0; k<NZ; ++k)
    {
      real_t x = (dx[0] * ((real_t)i + 0.5)) - L[0] / 2.0 ;
      real_t y = (dx[1] * ((real_t)j + 0.5)) - L[1] / 2.0 ;
      real_t z = (dx[2] * ((real_t)k + 0.5)) - L[2] / 2.0 ;

      real_t r = sqrt(pw2(x ) + pw2(y ) + pw2(z ));
      real_t phi = atan(y / x);
      real_t theta = acos(z / sqrt(pw2(x) + pw2(y) + pw2(z))); 

      real_t W = 0;
      real_t lap_Wr = 0, dx_K = 0, dy_K = 0, dz_K = 0;
      
      if(r < l)
        W = 0;
      else if(r<l+sigma)
      {
        W = std::pow(pow((r-l-sigma)/(sigma),6) - 1.0, 6 );
        // calculating derivative of W(r)/r
        lap_Wr = 90.0 * pow( (l-r+sigma) / sigma, 4.0)
          * pow(pow((r-l-sigma)/(sigma), 6) - 1.0, 4.0)
          * (-1.0+ 7 * pow( (l-r+sigma) / sigma, 6.0)) / (r * pw2(sigma));
        // calculating derivative of K
        dx_K = -36.0 * K_c * x * pow((l-r+sigma)/sigma,5.0)
          * pow(pow((r-l-sigma)/(sigma), 6) - 1.0 , 5) / (r * sigma);
        dy_K = -36.0 * K_c * y * pow((l-r+sigma)/sigma,5.0)
          * pow(pow((r-l-sigma)/(sigma), 6) - 1.0 , 5) / (r * sigma);
        dz_K = -36.0 * K_c * z * pow((l-r+sigma)/sigma,5.0)
          * pow(pow((r-l-sigma)/(sigma), 6) - 1.0 , 5) / (r * sigma);

      }
      else
        W = 1;
      real_t K = K_c * W;


      // initial guess
      X[0][INDEX(i,j,k)] = 1.0;

      X[1][INDEX(i,j,k)] = 0;
      X[2][INDEX(i,j,k)] = 0;
      X[3][INDEX(i,j,k)] = 0;
      
      // set coeficient of \Psi^5
      multigrid.setPolySrcAtPt(0, 16, i, j, k, -K*K / 12.0); 
      multigrid.setPolySrcAtPt(0, 17, i, j, k, -M * lap_Wr); 
          
      multigrid.setPolySrcAtPt(1, 4, i, j, k, -2.0 * dx_K / 3.0);
      multigrid.setPolySrcAtPt(2, 4, i, j, k, -2.0 * dy_K / 3.0);
      multigrid.setPolySrcAtPt(3, 4, i, j, k, -2.0 * dz_K / 3.0);

      multigrid.setShiftSrcAtPt(0, i, j, k, M / (2.0 * r) * (1 - W));

    }

    
    bd_handler->fillBoundary(X[0]._array, X[0].nx, X[0].ny, X[0].nz);
    bd_handler->fillBoundary(X[1]._array, X[1].nx, X[1].ny, X[1].nz);
    bd_handler->fillBoundary(X[2]._array, X[2].nx, X[2].ny, X[2].nz);
    bd_handler->fillBoundary(X[3]._array, X[3].nx, X[3].ny, X[3].nz);
    
    multigrid.initializeRhoHeirarchy();

      
    std::shared_ptr<tbox::HDFDatabase > hdf (new tbox::HDFDatabase("hdf_db"));

    std::string filename = "h5_init_data";

    mpi.Barrier();
    std::ifstream file(filename);
    if(file)
    {
      int rank = 0;
      while(rank < mpi.getSize())
      {
        if(rank == mpi.getRank())
        {
          hdf->open(filename, 1);
          const std::vector<double> & temp_X0 = hdf->getDoubleVector("X0");
          const std::vector<double> & temp_X1 = hdf->getDoubleVector("X1");
          const std::vector<double> & temp_X2 = hdf->getDoubleVector("X2");
          const std::vector<double> & temp_X3 = hdf->getDoubleVector("X3");

          // if file exist but corresponding database not exist
          if(temp_X0.empty())
            TBOX_ERROR("Getting empty array from file "<<filename<<"\n");

          tbox::pout<<"Read initial configuration database "<<"\n";
    
          for(int i = 0; i < temp_X0.size(); i++)
          {
            X[0]._array[i] = temp_X0[i];
            X[1]._array[i] = temp_X1[i];
            X[2]._array[i] = temp_X2[i];
            X[3]._array[i] = temp_X3[i];
          }
          flag = true;
          hdf->close();
        }
        mpi.Barrier();
        rank ++;
      }
    }
    else
    {
      multigrid.VCycles(num_vcycles);
      if(mpi.getRank() == 0)
      {
        // create and open the file
        hdf->create(filename);
        hdf->open(filename, 1);

        hdf->putDoubleArray("X0", X[0]._array, (NX+2*STENCIL_ORDER)*(NY+2*STENCIL_ORDER)*(NZ+2*STENCIL_ORDER));
        hdf->putDoubleArray("X1", X[1]._array, (NX+2*STENCIL_ORDER)*(NY+2*STENCIL_ORDER)*(NZ+2*STENCIL_ORDER));
        hdf->putDoubleArray("X2", X[2]._array, (NX+2*STENCIL_ORDER)*(NY+2*STENCIL_ORDER)*(NZ+2*STENCIL_ORDER));
        hdf->putDoubleArray("X3", X[3]._array, (NX+2*STENCIL_ORDER)*(NY+2*STENCIL_ORDER)*(NZ+2*STENCIL_ORDER));
        hdf->close();
      }
      flag = true;
      
    }

    // for(int i = 32; i < 64; i++)
    //   tbox::pout<<X[0][INDEX(i, 32, 32)]<<" ";
    // tbox::pout<<"\n\n\n";

    // for(int i = 32; i < 64; i++)
    //   tbox::pout<<X[1][INDEX(i, 32, 32)]<<" "<<"\n";
    // tbox::pout<<"\n";


    /**********start converting fields back*************/
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

    std::shared_ptr<pdat::CellData<real_t> > chi_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFchi_a_idx)));
    

    const hier::Box& box = chi_a_pdata->getGhostBox();
    


    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma12_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma12_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma13_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma13_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma22_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma23_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma23_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma33_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma33_a_idx)));


    std::shared_ptr<pdat::CellData<real_t> > A11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A12_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A12_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A13_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A13_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A22_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A23_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A23_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A33_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A33_a_idx)));

    std::shared_ptr<pdat::CellData<real_t> > DIFFK_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFK_a_idx)));

    arr_t DIFFchi_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        chi_a_pdata->getArrayData());

    arr_t DIFFgamma11_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma11_a_pdata->getArrayData());
    arr_t DIFFgamma12_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma12_a_pdata->getArrayData());
    arr_t DIFFgamma13_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma13_a_pdata->getArrayData());
    arr_t DIFFgamma22_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma22_a_pdata->getArrayData());
    arr_t DIFFgamma23_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma23_a_pdata->getArrayData());
    arr_t DIFFgamma33_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma33_a_pdata->getArrayData());

    arr_t A11_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A11_a_pdata->getArrayData());
    arr_t A12_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A12_a_pdata->getArrayData());
    arr_t A13_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A13_a_pdata->getArrayData());
    arr_t A22_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A22_a_pdata->getArrayData());
    arr_t A23_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A23_a_pdata->getArrayData());
    arr_t A33_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A33_a_pdata->getArrayData());

    arr_t DIFFK_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFK_a_pdata->getArrayData());
    
    const hier::Box& inner_box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    const int * inner_lower = &inner_box.lower()[0];
    const int * inner_upper = &inner_box.upper()[0];

    
    for(int k = inner_lower[2]; k <= inner_upper[2]; k++)
    {
      for(int j = inner_lower[1]; j <= inner_upper[1]; j++)
      {
        for(int i = inner_lower[0]; i <= inner_upper[0]; i++)
        {

          real_t x = (dx[0] * ((real_t)i + 0.5)) - L[0] / 2.0 ;
          real_t y = (dx[1] * ((real_t)j + 0.5)) - L[1] / 2.0 ;
          real_t z = (dx[2] * ((real_t)k + 0.5)) - L[2] / 2.0 ;

          real_t r = sqrt(pw2(x ) + pw2(y ) + pw2(z ));
          real_t phi = atan(y / x);
          real_t theta = acos(z / sqrt(pw2(x) + pw2(y) + pw2(z))); 

          real_t W = 0;
          if(r < l)
            W = 0;
          else if(r<l+sigma)
            W = std::pow(pow((r-l-sigma)/(sigma),6) - 1.0, 6 );
          else
            W = 1;
          real_t K = K_c * W;

          DIFFK_a(i, j, k) = K;
          DIFFchi_a(i,j,k) = 1.0 / pw2(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W)) - 1.0;

          real_t temp = 0.0;

          for(idx_t kk = 1; kk <=3; kk++)
            temp += multigrid.derivative(i, j, k, NX, NY, NZ, kk, X[kk]);
          
          A11_a(i, j, k) = std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) * (multigrid.derivative(i, j, k, NX, NY, NZ, 1, X[1]) + multigrid.derivative(i, j, k, NX, NY, NZ, 1, X[1]) - 2.0 * temp / 3.0);

          A12_a(i, j, k) = std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) * (multigrid.derivative(i, j, k, NX, NY, NZ, 1, X[2]) + multigrid.derivative(i, j, k, NX, NY, NZ, 2, X[1]));

          A13_a(i, j, k) = std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) * (multigrid.derivative(i, j, k, NX, NY, NZ, 1, X[3]) + multigrid.derivative(i, j, k, NX, NY, NZ, 3, X[1]));

          A22_a(i, j, k) = std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) * (multigrid.derivative(i, j, k, NX, NY, NZ, 2, X[2]) + multigrid.derivative(i, j, k, NX, NY, NZ, 2, X[2]) - 2.0 * temp / 3.0);
    
          A23_a(i, j, k) = std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) * (multigrid.derivative(i, j, k, NX, NY, NZ, 2, X[3]) + multigrid.derivative(i, j, k, NX, NY, NZ, 3, X[2]));

          A33_a(i, j, k) = std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) * (multigrid.derivative(i, j, k, NX, NY, NZ, 3, X[3]) + multigrid.derivative(i, j, k, NX, NY, NZ, 3, X[3]) - 2.0 * temp / 3.0);

        }
      }
    }
    // tbox::pout<<"\n";
    // for(int i = 32; i < 64; i++)
    //   tbox::pout<<A12_a(i, 32, 32)<<" "<<"\n";
    // tbox::pout<<"\n";

  }


}


// constructing kerr blackhole lattice initial data
// with conformal transverse-traceless decomposition
void bssn_ic_kerr_BHL_CTT(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln,
  real_t M,
  real_t a,
  real_t K_c,
  real_t relaxation_tolerance,
  idx_t num_vcycles,
  idx_t max_depth,
  double l,
  double sigma)
{
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  idx_t DIFFchi_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFchi"), variable_db->getContext("ACTIVE"));
  idx_t DIFFalpha_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFalpha"), variable_db->getContext("ACTIVE"));

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

  idx_t DIFFK_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFK"), variable_db->getContext("ACTIVE"));

  
  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
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


   for(int i = 0 ; i < 3; i++)
     L[i] = domain_upper[i] - domain_lower[i];


   std::string boundary_type = "periodic";
   multigridBdHandler * bd_handler = new multigridBdHandler(boundary_type, L, 10);

    
   idx_t NX = round(L[0] / dx[0]);
   idx_t NY = round(L[1] / dx[1]); 
   idx_t NZ = round(L[2] / dx[2]);


    // initializing solution vector
   CosmoArray<idx_t, real_t> * X = new CosmoArray<idx_t, real_t> [4];

   bool flag = false;

    X[0].init(NX, NY, NZ);
    X[1].init(NX, NY, NZ);
    X[2].init(NX, NY, NZ);
    X[3].init(NX, NY, NZ);

    idx_t molecule_n[4] = {27, 6, 6, 6};

    atom atom_tmp = {0};
    
    FASMultigrid multigrid(
      X, 4, molecule_n, max_depth, 2, relaxation_tolerance, L, NX, NY, NZ, bd_handler);

    /***********start initializaing equations **********************/

    multigrid.eqns[0][0].init(1, 1);
    //adding laplacian term
    atom_tmp.type =  multigrid.atom_type::lap;
    atom_tmp.u_id = 0;
    multigrid.eqns[0][0].add_atom(atom_tmp);

    //adding \Psi^-7 term
    for(int i = 1; i <= 3; i++)
      for(int j = 1; j <= 3; j++)
      {
        atom_tmp.type = i+1;
        atom_tmp.u_id = j;
        if(i == j)
          multigrid.eqns[0][3*(i-1)+j].init(3, 0.5 - 1.0/6.0);
        else
          multigrid.eqns[0][3*(i-1)+j].init(3, 0.25);
        multigrid.eqns[0][3*(i-1)+j].add_atom(atom_tmp);

        multigrid.eqns[0][3*(i-1)+j].add_atom(atom_tmp);

        atom_tmp.type = 1;
        atom_tmp.u_id = 0;
        atom_tmp.value = -7;
        multigrid.eqns[0][3*(i-1)+j].add_atom(atom_tmp);
      }
    for(int i = 1; i <= 3; i++)
      for(int j = i+1; j <= 3; j++)
      {
        atom_tmp.type = i+1;
        atom_tmp.u_id = j;
        multigrid.eqns[0][(9 + i+j - 2)].init(3, 0.5);
        multigrid.eqns[0][(9 + i+j - 2)].add_atom(atom_tmp);
      
        atom_tmp.type = j+1;
        atom_tmp.u_id = i;
        multigrid.eqns[0][(9 + i+j - 2)].add_atom(atom_tmp);

        atom_tmp.type = 1;
        atom_tmp.u_id = 0;
        atom_tmp.value = -7;
        multigrid.eqns[0][(9 + i+j - 2)].add_atom(atom_tmp);
      }

    for(int i = 1; i <= 3; i++)
      for(int j = i+1; j <= 3; j++)
      {
        multigrid.eqns[0][(12 + i+j - 2)].init(3, -1.0/3.0);
        atom_tmp.type = i+1;
        atom_tmp.u_id = i;
        multigrid.eqns[0][(12+i+j-2)].add_atom(atom_tmp);

        atom_tmp.type = j+1;
        atom_tmp.u_id = j;
        multigrid.eqns[0][(12+i+j-2)].add_atom(atom_tmp);

        atom_tmp.type = 1;
        atom_tmp.u_id = 0;
        atom_tmp.value = -7;
        multigrid.eqns[0][(12 + i+j - 2)].add_atom(atom_tmp);
      }

    // adding \partial X \partial S terms
    // where S is divergen shift
    for(int i = 1; i <= 3; i++)
      for(int j = 1; j <= 3; j++)
      {
        if( i == 3 && j == 3) continue;

        multigrid.eqns[0][15+3*(i-1)+j].init(2, 1.0);
        
        atom_tmp.type = i+1;
        atom_tmp.u_id = j;
        multigrid.eqns[0][15+3*(i-1)+j].add_atom(atom_tmp);

        atom_tmp.type = 1;
        atom_tmp.u_id = 0;
        atom_tmp.value = -7;
        multigrid.eqns[0][15+3*(i-1)+j].add_atom(atom_tmp);
      }

    // adding pure \partial S terms
    multigrid.eqns[0][24].init(1, 1.0);
    atom_tmp.type = 1;
    atom_tmp.u_id = 0;
    atom_tmp.value = -7;
    multigrid.eqns[0][24].add_atom(atom_tmp);

    multigrid.eqns[0][25].init(1, 1.0);

    // adding \Psi^5 terms
    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 5;
    multigrid.eqns[0][25].add_atom(atom_tmp);

    // adding Lap(divergence part) as a ^0 polynomial
    multigrid.eqns[0][26].init(1, 1.0);
    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 0;
    multigrid.eqns[0][26].add_atom(atom_tmp);
    
    // start adding momentum constraints
    multigrid.eqns[1][0].init(1, 1.0);
    multigrid.eqns[1][1].init(1, 1.0/3.0);
    multigrid.eqns[1][2].init(1, 1.0/3.0);
    multigrid.eqns[1][3].init(1, 1.0/3.0);
    
    multigrid.eqns[1][4].init(1, 1.0);
    multigrid.eqns[1][5].init(1, 1.0);

    multigrid.eqns[2][0].init(1, 1.0);
    multigrid.eqns[2][1].init(1, 1.0/3.0);
    multigrid.eqns[2][2].init(1, 1.0/3.0);
    multigrid.eqns[2][3].init(1, 1.0/3.0);
    multigrid.eqns[2][4].init(1, 1.0);
    multigrid.eqns[2][5].init(1, 1.0);


    multigrid.eqns[3][0].init(1, 1.0);
    multigrid.eqns[3][1].init(1, 1.0/3.0);
    multigrid.eqns[3][2].init(1, 1.0/3.0);
    multigrid.eqns[3][3].init(1, 1.0/3.0);
    multigrid.eqns[3][4].init(1, 1.0);
    multigrid.eqns[3][5].init(1, 1.0);



    //adding terms to eqn 1
    atom_tmp.type = multigrid.atom_type::lap;
    atom_tmp.u_id = 1;
    multigrid.eqns[1][0].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der11;
    atom_tmp.u_id = 1;
    multigrid.eqns[1][1].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der12;
    atom_tmp.u_id = 2;
    multigrid.eqns[1][2].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der13;
    atom_tmp.u_id = 3;
    multigrid.eqns[1][3].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 6;
    multigrid.eqns[1][4].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 0;
    multigrid.eqns[1][5].add_atom(atom_tmp);

  
    //adding terms to eqn 2
    atom_tmp.type = multigrid.atom_type::lap;
    atom_tmp.u_id = 2;
    multigrid.eqns[2][0].add_atom(atom_tmp);


    atom_tmp.type = multigrid.atom_type::der12;
    atom_tmp.u_id = 1;
    multigrid.eqns[2][1].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der22;
    atom_tmp.u_id = 2;
    multigrid.eqns[2][2].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der23;
    atom_tmp.u_id = 3;
    multigrid.eqns[2][3].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 6;
    multigrid.eqns[2][4].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 0;
    multigrid.eqns[2][5].add_atom(atom_tmp);

    
    //adding terms to eqn 3
    atom_tmp.type = multigrid.atom_type::lap;
    atom_tmp.u_id = 3;
    multigrid.eqns[3][0].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der13;
    atom_tmp.u_id = 1;
    multigrid.eqns[3][1].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der23;
    atom_tmp.u_id = 2;
    multigrid.eqns[3][2].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::der33;
    atom_tmp.u_id = 3;
    multigrid.eqns[3][3].add_atom(atom_tmp);

    
    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 0;
    multigrid.eqns[3][5].add_atom(atom_tmp);

    atom_tmp.type = multigrid.atom_type::poly;
    atom_tmp.u_id = 0;
    atom_tmp.value = 6;
    multigrid.eqns[3][4].add_atom(atom_tmp);

  
  /*********ending initializing equations************/
    
    for(int i=0; i<NX; ++i) 
    for(int j=0; j<NY; ++j) 
    for(int k=0; k<NZ; ++k)
    {
      real_t x = (dx[0] * ((real_t)i + 0.5)) - L[0] / 2.0 ;
      real_t y = (dx[1] * ((real_t)j + 0.5)) - L[1] / 2.0 ;
      real_t z = (dx[2] * ((real_t)k + 0.5)) - L[2] / 2.0 ;

      real_t r = sqrt(pw2(x ) + pw2(y ) + pw2(z ));
      real_t phi = atan(y / x);
      real_t theta = acos(z / sqrt(pw2(x) + pw2(y) + pw2(z))); 

      real_t W = 0;
      real_t lap_Wr = 0, dx_K = 0, dy_K = 0, dz_K = 0;
      real_t sup_X1 = 0, sup_X2 = 0, sup_X3 = 0;

      real_t ds = (l-r+sigma)/sigma;

      real_t d1S1 = 0, d1S2 = 0, d1S3 = 0, d2S1 = 0, d2S2 = 0,
        d2S3 = 0, d3S1 = 0, d3S2 = 0, d3S3 = 0;
      
      if(r < l)
      {
        W = 0;
        d1S1 = - 3.0 * a * x * y / pow(r, 5.0);
        d1S2 = - a * (r*r - 3 * x * x) / pow(r, 5.0);
        d2S1 = a * (r*r - 3 * y * y) / pow(r, 5.0);
        d2S2 = 3.0 * a * x * y / pow(r, 5.0);
        d3S1 = -3.0 * a * z * y / pow(r, 5.0);
        d3S2 = 3.0 * a * x * z / pow(r, 5.0);
      }
      else if(r<l+sigma)
      {
        W = std::pow(pow((r-l-sigma)/(sigma),6) - 1.0, 6 );
        // calculating derivative of W(r)/r
        lap_Wr = 90.0 * pow( (l-r+sigma) / sigma, 4.0)
          * pow(pow((r-l-sigma)/(sigma), 6) - 1.0, 4.0)
          * (-1.0+ 7 * pow( (l-r+sigma) / sigma, 6.0)) / (r * pw2(sigma));
        // calculating derivative of K
        dx_K = -36.0 * K_c * x * pow((l-r+sigma)/sigma,5.0)
          * pow(pow((r-l-sigma)/(sigma), 6) - 1.0 , 5) / (r * sigma);
        dy_K = -36.0 * K_c * y * pow((l-r+sigma)/sigma,5.0)
          * pow(pow((r-l-sigma)/(sigma), 6) - 1.0 , 5) / (r * sigma);
        dz_K = -36.0 * K_c * z * pow((l-r+sigma)/sigma,5.0)
          * pow(pow((r-l-sigma)/(sigma), 6) - 1.0 , 5) / (r * sigma);

        sup_X1 = - 36.0 * a * y * pow(ds, 4) * pow( 1.0 - pow(ds, 6), 4)
          * (30.0 * r * pow(ds, 6) - 5.0 * r * (1- pow(ds, 6))
             - 2.0 * (l-r+sigma) * (1.0 - pow(ds, 6.0))) / pow(r, 4.0) / sigma / sigma; 

        sup_X2 = 36.0 * a * x * pow(ds, 4) * pow( 1.0 - pow(ds, 6), 4)
          * (30.0 * r * pow(ds, 6) - 5.0 * r * (1- pow(ds, 6))
             - 2.0 * (l-r+sigma) * (1.0 - pow(ds, 6.0))) / pow(r, 4.0) / sigma / sigma;
        
        sup_X3 = 0;


        d1S1 = 3*a*x*y*pow(r,-5)*(-1 - 12*r*pow(sigma,-6)*pow(l - r + sigma,5)*pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),5) + pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),6));
        d1S2 = a*pow(r,-5)*(36*r*pow(sigma,-6)*pow(l - r + sigma,5)*pow(x,2)*pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),5) - pow(r,2)*(1 - pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),6)) + 
                            3*pow(x,2)*(1 - pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),6)));
        d2S1 = -(a*pow(r,-5)*(36*r*pow(sigma,-6)*pow(l - r + sigma,5)*pow(y,2)*pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),5) - pow(r,2)*(1 - pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),6)) + 
                              3*pow(y,2)*(1 - pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),6))));
        d2S2 = -3*a*x*y*pow(r,-5)*(-1 - 12*r*pow(sigma,-6)*pow(l - r + sigma,5)*pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),5) + pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),6));

        d3S1 = 3*a*y*z*pow(r,-5)*(-1 - 12*r*pow(sigma,-6)*pow(l - r + sigma,5)*pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),5) + pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),6));

        d3S2 = -3*a*x*z*pow(r,-5)*pow(sigma,-36)*(-pow(sigma,36) - 12*r*pow(l - r + sigma,5)*pow(pow(sigma,6) - pow(l - r + sigma,6),5) + pow(pow(sigma,6) - pow(l - r + sigma,6),6));
        
      }
      else
        W = 1;
      

      real_t K = K_c * W;


      // initial guess
      X[0][INDEX(i,j,k)] = 1.0;

      X[1][INDEX(i,j,k)] = 0;
      X[2][INDEX(i,j,k)] = 0;
      X[3][INDEX(i,j,k)] = 0;
      
      // set coeficient of \Psi^5
      multigrid.setPolySrcAtPt(0, 25, i, j, k, -K*K / 12.0);
      
      multigrid.setPolySrcAtPt(0, 26, i, j, k, -M * lap_Wr);
      real_t const_f_of_psi7;
      const_f_of_psi7 = 
        + (2.0 * (d1S1*d1S1 + d2S2*d2S2 + d3S3*d3S3)
        + d1S2*d1S2 + d2S1*d2S1 + d1S3*d1S3 + d3S1*d3S1 + d2S3*d2S3 + d3S2*d3S2
        + d1S2*d2S1 + d2S1*d1S2 + d1S3*d3S1 + d3S1*d1S3 + d2S3*d3S2 + d3S2*d2S3)/4.0;
      
      multigrid.setPolySrcAtPt(0, 24, i, j, k, const_f_of_psi7); 

      multigrid.setPolySrcAtPt(0, 16, i, j, k, d1S1); // coef of d1X1
      multigrid.setPolySrcAtPt(0, 17, i, j, k, 0.5*(d1S2 + d2S1)); // coef of d1X2
      multigrid.setPolySrcAtPt(0, 18, i, j, k, 0.5*d3S1); // coef of d1X3

      multigrid.setPolySrcAtPt(0, 19, i, j, k, 0.5*(d2S1 + d1S2)); // coef of d2X1
      multigrid.setPolySrcAtPt(0, 20, i, j, k, d2S2); // coef of d2X2
      multigrid.setPolySrcAtPt(0, 21, i, j, k, 0.5*d3S2); // coef of d2X3

      multigrid.setPolySrcAtPt(0, 22, i, j, k, 0.5*d3S1); // coef of d3X1
      multigrid.setPolySrcAtPt(0, 23, i, j, k, 0.5*d3S2); // coef of d3X2
      
      multigrid.setPolySrcAtPt(1, 4, i, j, k, -2.0 * dx_K / 3.0);
      multigrid.setPolySrcAtPt(2, 4, i, j, k, -2.0 * dy_K / 3.0); 
      multigrid.setPolySrcAtPt(3, 4, i, j, k, -2.0 * dz_K / 3.0);

      multigrid.setPolySrcAtPt(1, 5, i, j, k, sup_X1);
      multigrid.setPolySrcAtPt(2, 5, i, j, k, sup_X2);
      multigrid.setPolySrcAtPt(3, 5, i, j, k, sup_X3);
      
      multigrid.setShiftSrcAtPt(0, i, j, k, M / (2.0 * r) * (1 - W));
      
      multigrid.setShiftSrcAtPt(1, i, j, k, a * y * (1 - W) / pw3(r));
      multigrid.setShiftSrcAtPt(2, i, j, k, - a * x * (1 - W) / pw3(r));
      multigrid.setShiftSrcAtPt(3, i, j, k, 0);
    }

    bd_handler->fillBoundary(X[0]._array, X[0].nx, X[0].ny, X[0].nz);
    bd_handler->fillBoundary(X[1]._array, X[1].nx, X[1].ny, X[1].nz);
    bd_handler->fillBoundary(X[2]._array, X[2].nx, X[2].ny, X[2].nz);
    bd_handler->fillBoundary(X[3]._array, X[3].nx, X[3].ny, X[3].nz);
    
    multigrid.initializeRhoHeirarchy();

      
    std::shared_ptr<tbox::HDFDatabase > hdf (new tbox::HDFDatabase("hdf_db"));

    std::string filename = "h5_init_data";

    mpi.Barrier();
    std::ifstream file(filename);
    
    // if file exists, reading from file
    // run multigrid solver otherwise
    if(file)
    {
      int rank = 0;
      while(rank < mpi.getSize())
      {
        if(rank == mpi.getRank())
        {
          hdf->open(filename, 1);
          const std::vector<double> & temp_X0 = hdf->getDoubleVector("X0");
          const std::vector<double> & temp_X1 = hdf->getDoubleVector("X1");
          const std::vector<double> & temp_X2 = hdf->getDoubleVector("X2");
          const std::vector<double> & temp_X3 = hdf->getDoubleVector("X3");

          // if file exist but corresponding database not exist
          if(temp_X0.empty())
            TBOX_ERROR("Getting empty array from file "<<filename<<"\n");

          tbox::pout<<"Read initial configuration database "<<"\n";
    
          for(int i = 0; i < temp_X0.size(); i++)
          {
            X[0]._array[i] = temp_X0[i];
            X[1]._array[i] = temp_X1[i];
            X[2]._array[i] = temp_X2[i];
            X[3]._array[i] = temp_X3[i];
          }
          flag = true;
          hdf->close();
        }
        mpi.Barrier();
        rank ++;
      }
      //      multigrid.VCycles(num_vcycles);
    }
    else
    {
      multigrid.VCycles(num_vcycles);
      if(mpi.getRank() == 0)
      {
        // create and open the file
        hdf->create(filename);
        hdf->open(filename, 1);

        hdf->putDoubleArray("X0", X[0]._array, (NX+2*STENCIL_ORDER)*(NY+2*STENCIL_ORDER)*(NZ+2*STENCIL_ORDER));
        hdf->putDoubleArray("X1", X[1]._array, (NX+2*STENCIL_ORDER)*(NY+2*STENCIL_ORDER)*(NZ+2*STENCIL_ORDER));
        hdf->putDoubleArray("X2", X[2]._array, (NX+2*STENCIL_ORDER)*(NY+2*STENCIL_ORDER)*(NZ+2*STENCIL_ORDER));
        hdf->putDoubleArray("X3", X[3]._array, (NX+2*STENCIL_ORDER)*(NY+2*STENCIL_ORDER)*(NZ+2*STENCIL_ORDER));
        hdf->close();
      }
      flag = true;
      
    }
    /**************debuging *************/
    

    // for(int i = NX/2; i < NX; i++)
    //   tbox::pout<<(X[1][INDEX(i, NX/2, NX/2)] + X[1][INDEX(i, NX/2-1, NX/2)] +
    //                X[1][INDEX(i, NX/2, NX/2-1)] + X[1][INDEX(i, NX/2-1, NX/2-1)]) / 4.0<<", ";

    // tbox::pout<<"\n\n\n";

    //    tbox::pout<<"\n\n\n";


    
    // for(int i = NX/2; i <NX; i++)
    //   tbox::pout<<(X[3][INDEX(i, NX/2, NX/2)] + X[3][INDEX(i, NX/2-1, NX/2)] +
    //                X[3][INDEX(i, NX/2, NX/2-1)] + X[3][INDEX(i, NX/2-1, NX/2-1)]) / 4.0<<", ";
    // tbox::pout<<"\n\n\n";

    /*********ending debuging part*************/

    /**********start converting fields back to BSSN fields*************/
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

    std::shared_ptr<pdat::CellData<real_t> > chi_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFchi_a_idx)));

    std::shared_ptr<pdat::CellData<real_t> > DIFFalpha_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFalpha_a_idx)));


    const hier::Box& box = chi_a_pdata->getGhostBox();
    


    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma12_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma12_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma13_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma13_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma22_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma23_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma23_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma33_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma33_a_idx)));


    std::shared_ptr<pdat::CellData<real_t> > A11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A12_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A12_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A13_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A13_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A22_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A23_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A23_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A33_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A33_a_idx)));

    std::shared_ptr<pdat::CellData<real_t> > DIFFK_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFK_a_idx)));

    arr_t DIFFchi_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        chi_a_pdata->getArrayData());
    arr_t DIFFalpha_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFalpha_a_pdata->getArrayData());

    
    arr_t DIFFgamma11_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma11_a_pdata->getArrayData());
    arr_t DIFFgamma12_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma12_a_pdata->getArrayData());
    arr_t DIFFgamma13_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma13_a_pdata->getArrayData());
    arr_t DIFFgamma22_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma22_a_pdata->getArrayData());
    arr_t DIFFgamma23_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma23_a_pdata->getArrayData());
    arr_t DIFFgamma33_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma33_a_pdata->getArrayData());

    arr_t A11_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A11_a_pdata->getArrayData());
    arr_t A12_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A12_a_pdata->getArrayData());
    arr_t A13_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A13_a_pdata->getArrayData());
    arr_t A22_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A22_a_pdata->getArrayData());
    arr_t A23_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A23_a_pdata->getArrayData());
    arr_t A33_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A33_a_pdata->getArrayData());

    arr_t DIFFK_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFK_a_pdata->getArrayData());
    
    const hier::Box& inner_box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    const int * inner_lower = &inner_box.lower()[0];
    const int * inner_upper = &inner_box.upper()[0];

    
    for(int k = inner_lower[2]; k <= inner_upper[2]; k++)
    {
      for(int j = inner_lower[1]; j <= inner_upper[1]; j++)
      {
        for(int i = inner_lower[0]; i <= inner_upper[0]; i++)
        {

          real_t x = (dx[0] * ((real_t)i + 0.5)) - L[0] / 2.0 ;
          real_t y = (dx[1] * ((real_t)j + 0.5)) - L[1] / 2.0 ;
          real_t z = (dx[2] * ((real_t)k + 0.5)) - L[2] / 2.0 ;

          real_t r = sqrt(pw2(x ) + pw2(y ) + pw2(z ));
          real_t phi = atan(y / x);
          real_t theta = acos(z / sqrt(pw2(x) + pw2(y) + pw2(z))); 

          real_t W = 0;
          real_t sA11 = 0, sA12 = 0, sA13 = 0, sA22 = 0, sA33 = 0, sA23 = 0;

          if(r < l)
          {
            W = 0;
            sA11 = -6.0 * a * x * y/ pow(r, 5.0);
            sA12 = 3.0 * a * (pw2(x) - pw2(y)) / pow(r, 5.0);
            sA13 = -3.0 * a * y * z/ pow(r, 5.0);
            sA22 = 6.0 * a * x * y/ pow(r, 5.0);
            sA23 = 3.0 * a * x * z/ pow(r, 5.0);
          }
          else if(r<l+sigma)
          {
            W = std::pow(pow((r-l-sigma)/(sigma),6) - 1.0, 6 );
            sA11 = 
              +6*a*x*y*pow(r,-5)*(-1 - 12*r*pow(sigma,-6)*pow(l - r + sigma,5)*pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),5) + pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),6));
            sA12 =                
               + 3*a*pow(r,-5)*(12*r*pow(sigma,-6)*pow(l - r + sigma,5)*pow(x,2)*pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),5) - 12*r*pow(sigma,-6)*pow(l - r + sigma,5)*pow(y,2)*pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),5) + 
                                pow(x,2)*(1 - pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),6)) - pow(y,2)*(1 - pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),6)));
            sA13 = +3*a*y*z*pow(r,-5)*(-1 - 12*r*pow(sigma,-6)*pow(l - r + sigma,5)*pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),5) + pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),6));
            sA22 =              
              -6*a*x*y*pow(r,-5)*(-1 - 12*r*pow(sigma,-6)*pow(l - r + sigma,5)*pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),5) + pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),6));
            sA23 = -3*a*x*z*pow(r,-5)*(-1 - 12*r*pow(sigma,-6)*pow(l - r + sigma,5)*pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),5) + pow(1 - pow(sigma,-6)*pow(l - r + sigma,6),6));
          }
          else
            W = 1;
          real_t K = K_c * W;

          real_t ds = (l-r+sigma)/sigma;
          
          DIFFK_a(i, j, k) = K;
          DIFFchi_a(i,j,k) = 1.0 / pw2(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W)) - 1.0;
          //          DIFFalpha_a(i, j, k) = DIFFchi_a(i, j, k);
          
          
          real_t temp = 0.0;

          for(idx_t kk = 1; kk <=3; kk++)
            temp += multigrid.derivative(i, j, k, NX, NY, NZ, kk, X[kk]);

          
          A11_a(i, j, k) = std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0)
            * (multigrid.derivative(i, j, k, NX, NY, NZ, 1, X[1])
               + multigrid.derivative(i, j, k, NX, NY, NZ, 1, X[1])
               - 2.0 * temp / 3.0);

          // if(i == (inner_upper[1] + 1)/2 && j == (inner_upper[1] + 1)/2 && k ==(inner_upper[1] + 1)/2)
          //   std::cout<<"!@!!! "<<X[0][INDEX(i, j, k)]<<"\n";
          
          A12_a(i, j, k) = std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0)
            * (multigrid.derivative(i, j, k, NX, NY, NZ, 1, X[2])
               + multigrid.derivative(i, j, k, NX, NY, NZ, 2, X[1]));

          
          A13_a(i, j, k) = std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) *
            (multigrid.derivative(i, j, k, NX, NY, NZ, 1, X[3])
             + multigrid.derivative(i, j, k, NX, NY, NZ, 3, X[1]));

          
          A22_a(i, j, k) = std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) *
            (multigrid.derivative(i, j, k, NX, NY, NZ, 2, X[2])
             + multigrid.derivative(i, j, k, NX, NY, NZ, 2, X[2])
             - 2.0 * temp / 3.0);
          
          A23_a(i, j, k) = std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) *
            (multigrid.derivative(i, j, k, NX, NY, NZ, 2, X[3])
             + multigrid.derivative(i, j, k, NX, NY, NZ, 3, X[2]));
          
          A33_a(i, j, k) = std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) *
            (multigrid.derivative(i, j, k, NX, NY, NZ, 3, X[3])
             + multigrid.derivative(i, j, k, NX, NY, NZ, 3, X[3])
             - 2.0 * temp / 3.0);

          A11_a(i, j, k) += std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) * sA11;
          A12_a(i, j, k) += std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) * sA12;
          A13_a(i, j, k) += std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) * sA13;
          A22_a(i, j, k) += std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) * sA22;
          A23_a(i, j, k) += std::pow(X[0][INDEX(i, j, k)] + M / (2.0 * r) * (1 - W), -6.0) * sA23;
        }
      }
    }

    
    // for(int i = (inner_upper[0]+1)/2; i <= inner_upper[0]; i++)
    // {
    //   int j = (inner_upper[1] + 1) / 2, k = (inner_upper[2] +1)/ 2;
    //   std::cout<<(DIFFchi_a(i, j, k) + DIFFchi_a(i, j-1, k)
    //               + DIFFchi_a(i, j, k-1) + DIFFchi_a(i, j-1, k-1) ) /4.0
    //            <<",";
    // }
    // std::cout<<"\n\n\n";

    
    // for(int i = (inner_upper[0]+1)/2; i <= inner_upper[0]; i++)
    // {
    //   int j = (inner_upper[1] + 1) / 2, k = (inner_upper[2] +1)/ 2;
    //   std::cout<<(A11_a(i, j, k) + A11_a(i, j-1, k)
    //               + A11_a(i, j, k-1) + A11_a(i, j-1, k-1) ) /4.0
    //            <<",";
    // }
    // std::cout<<"\n\n\n";

    // for(int i = (inner_upper[0]+1)/2; i <= inner_upper[0]; i++)
    // {
    //   int j = (inner_upper[1] + 1) / 2, k = (inner_upper[2] +1)/ 2;
    //   std::cout<<(A22_a(i, j, k) + A22_a(i, j-1, k)
    //               + A22_a(i, j, k-1) + A22_a(i, j-1, k-1) ) /4.0
    //            <<",";
    // }
    // std::cout<<"\n\n\n";


    // for(int i = (inner_upper[0]+1)/2; i <= inner_upper[0]; i++)
    // {
    //   int j = (inner_upper[1] + 1) / 2, k = (inner_upper[2] +1)/ 2;
    //   std::cout<<(A33_a(i, j, k) + A33_a(i, j-1, k)
    //               + A33_a(i, j, k-1) + A33_a(i, j-1, k-1) ) /4.0
    //            <<",";
    // }
    // std::cout<<"\n\n\n";

    // for(int i = (inner_upper[0]+1)/2; i <= inner_upper[0]; i++)
    // {
    //   int j = (inner_upper[1] + 1) / 2, k = (inner_upper[2] +1)/ 2;
    //   std::cout<<(A12_a(i, j, k) + A12_a(i, j-1, k)
    //               + A12_a(i, j, k-1) + A12_a(i, j-1, k-1) ) /4.0
    //            <<",";
    // }
    // std::cout<<"\n\n\n";

    // for(int i = (inner_upper[0]+1)/2; i <= inner_upper[0]; i++)
    // {
    //   //      int j = (inner_upper[1] + 5) / 2, k = (inner_upper[2] +5)/ 2;
    //   std::cout<<(A13_a(i, i, i)) /4.0
    //            <<",";
    // }
    // std::cout<<"\n\n\n";

    // for(int i = (inner_upper[0]+1)/2; i <= inner_upper[0]; i++)
    // {
    //   //      int j = (inner_upper[1] + 5) / 2, k = (inner_upper[2] +5)/ 2;
    //   std::cout<<(A23_a(i, i, i)) /4.0
    //            <<",";
    // }
    // std::cout<<"\n\n\n";

    // for(int i = (inner_upper[0]+1)/2; i <= inner_upper[0]; i++)
    // {
    //   int j = (inner_upper[1] + 5) / 2, k = (inner_upper[2] +5)/ 2;
    //   std::cout<<(A23_a(i, j, k) + A23_a(i, j-1, k)
    //               + A23_a(i, j, k-1) + A23_a(i, j-1, k-1) ) /4.0
    //            <<",";
    // }
    // std::cout<<"\n\n\n";




  }


}

  //http://arxiv.org/abs/1001.4077v1
void bssn_ic_kerr_blackhole(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln)
{
# if ! USE_BSSN_SHIFT
  std::cerr << "Waning! USE_BSSN_SHIFT is suggested to be enabled in blackhole test " << std::endl;
# endif

  // spin of momentum
  real_t a = 0.6;

  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

  idx_t DIFFchi_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFchi"), variable_db->getContext("ACTIVE"));

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

  idx_t DIFFK_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFK"), variable_db->getContext("ACTIVE"));

  
  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
       hierarchy->getGridGeometry()));
   TBOX_ASSERT(grid_geometry_);
   geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;


    const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

    std::shared_ptr<pdat::CellData<real_t> > chi_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFchi_a_idx)));
    

    const hier::Box& box = chi_a_pdata->getGhostBox();
    


    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma12_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma12_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma13_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma13_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma22_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma23_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma23_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma33_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFgamma33_a_idx)));


    std::shared_ptr<pdat::CellData<real_t> > A11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A12_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A12_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A13_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A13_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A22_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A23_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A23_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A33_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(A33_a_idx)));

    std::shared_ptr<pdat::CellData<real_t> > DIFFK_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFK_a_idx)));

    arr_t chi_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        chi_a_pdata->getArrayData());

    arr_t DIFFgamma11_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma11_a_pdata->getArrayData());
    arr_t DIFFgamma12_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma12_a_pdata->getArrayData());
    arr_t DIFFgamma13_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma13_a_pdata->getArrayData());
    arr_t DIFFgamma22_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma22_a_pdata->getArrayData());
    arr_t DIFFgamma23_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma23_a_pdata->getArrayData());
    arr_t DIFFgamma33_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFgamma33_a_pdata->getArrayData());

    arr_t A11_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A11_a_pdata->getArrayData());
    arr_t A12_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A12_a_pdata->getArrayData());
    arr_t A13_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A13_a_pdata->getArrayData());
    arr_t A22_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A22_a_pdata->getArrayData());
    arr_t A23_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A23_a_pdata->getArrayData());
    arr_t A33_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        A33_a_pdata->getArrayData());

    arr_t DIFFK_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        DIFFK_a_pdata->getArrayData());

    
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

          real_t r = sqrt(pw2(x ) + pw2(y ) + pw2(z ));
          real_t phi = atan(y / x);
          real_t theta = acos(z / sqrt(pw2(x) + pw2(y) + pw2(z))); 
          
          
          real_t rp = 1 + sqrt(1 - pw2(a));
          real_t rm = 1 - sqrt(1 - pw2(a));
          real_t r_BL = r * pw2(1 + rp / (4.0 * r));
          real_t Sigma = pw2(r_BL) + pw2(a) * pw2(cos(theta));
          real_t Delta = pw2(r_BL) - 2.0 * r_BL + pw2(a);
          real_t A = pw2(pw2(r_BL) + pw2(a))
            - Delta * pw2(a * sin(theta));
          
          real_t ms11 = Sigma * pw2(r + rp / 4.0 ) / (pw3(r) * (r_BL - rm));
          real_t ms22 = Sigma;
          real_t ms33 = A * pw2(sin(theta)) / Sigma;

          real_t Ks13 = a * pw2(sin(theta)) / (Sigma * sqrt(Sigma  * A * r * (r_BL - rm)))
            * (1 + rp / (4.0 * r)) * (3.0 * pw2(pw2(r_BL)) + 2.0 * pw2(a) * pw2(r_BL)
                                      - pw2(pw2(a)) - pw2(a * sin(theta)) * (pw2(r_BL) - pw2(a)));
          real_t Ks23 = - (r - rp / (4.0 )) * sqrt( (r_BL - rm)/r ) * 2.0 * a * a * a
            * r_BL * cos(theta) * pw2(sin(theta)) * sin(theta) / (Sigma * sqrt(A * Sigma));


          
          DIFFgamma11_a(i, j, k) = (pw2(x)*ms11)/pw2(r) + 
            ((pw2(x)*(pw2(x) + pw2(y))*pw2(z)*ms22)/pow(r,4) 
             +  pw2(y)*ms33)/pow(pw2(x) + pw2(y),2);

          DIFFgamma12_a(i, j, k) = (x*y*(pw2(r)*pow(pw2(x) + pw2(y),2)*ms11 + 
       (pw2(x) + pw2(y))*pw2(z)*ms22 - pow(r,4)*ms33))/
            (pow(r,4)*pow(pw2(x) + pw2(y),2));

          DIFFgamma13_a(i, j, k) = (x*z*(pw2(r)*ms11 - ms22))/pow(r,4);

          DIFFgamma22_a(i, j, k)  = (pw2(y)*ms11)/pw2(r) + 
            ((pw2(y)*(pw2(x) + pw2(y))*pw2(z)*ms22)/pow(r,4) 
             + pw2(x)*ms33)/pow(pw2(x) + pw2(y),2);
          DIFFgamma23_a(i, j, k) = (y*z*(pw2(r)*ms11 - ms22))/pow(r,4);

          DIFFgamma33_a(i, j, k) = (pw2(r)*pw2(z)*ms11 + (pw2(x) + pw2(y))*ms22)/
   pow(r,4);
          
          real_t det = DIFFgamma11_a(i, j, k) * DIFFgamma22_a(i, j, k) * DIFFgamma33_a(i, j, k)
            + DIFFgamma12_a(i, j, k) * DIFFgamma23_a(i, j, k) * DIFFgamma13_a(i, j, k)
            + DIFFgamma12_a(i, j, k) * DIFFgamma23_a(i, j, k) * DIFFgamma13_a(i, j, k)
            - DIFFgamma13_a(i, j, k) * DIFFgamma22_a(i, j, k) * DIFFgamma13_a(i, j, k)
            - DIFFgamma12_a(i, j, k) * DIFFgamma12_a(i, j, k) * DIFFgamma33_a(i, j, k)
            - DIFFgamma23_a(i, j, k) * DIFFgamma23_a(i, j, k) * DIFFgamma11_a(i, j, k);

          chi_a(i,j,k) = pow(det, - 1/6.0);

          
          real_t gammai11 = (DIFFgamma22_a(i, j, k) * DIFFgamma33_a(i, j, k) - pw2(DIFFgamma23_a(i, j, k))) / det;
          real_t gammai22 = (DIFFgamma11_a(i, j, k) * DIFFgamma33_a(i, j, k) - pw2(DIFFgamma13_a(i, j, k))) / det;
          real_t gammai33 = (DIFFgamma11_a(i, j, k) * DIFFgamma22_a(i, j, k) - pw2(DIFFgamma12_a(i, j, k))) / det;
          real_t gammai12 = (DIFFgamma13_a(i, j, k) * DIFFgamma23_a(i, j, k) - DIFFgamma12_a(i, j, k)*(DIFFgamma33_a(i, j, k))) / det;
          real_t gammai13 = (DIFFgamma12_a(i, j, k) * DIFFgamma23_a(i, j, k) - DIFFgamma13_a(i, j, k)*(DIFFgamma22_a(i, j, k))) / det;
          real_t gammai23 = (DIFFgamma12_a(i, j, k) * DIFFgamma13_a(i, j, k) - DIFFgamma23_a(i, j, k)*(DIFFgamma11_a(i, j, k))) / det;


          real_t K11 = (-2*x*y*((pw2(x) + pw2(y))*Ks13 + 
                                sqrt((pw2(x) + pw2(y))/pw2(r))*z*Ks23))/
            (sqrt(pw2(r))*pow(pw2(x) + pw2(y),2));
          
          real_t K12 = ((x - y)*(x + y)*((pw2(x) + pw2(y))*Ks13 + 
                                         sqrt((pw2(x) + pw2(y))/pw2(r))*z*Ks23))/
            (sqrt(pw2(r))*pow(pw2(x) + pw2(y),2));
          
          real_t K13 = (y*(-(z*Ks13) + sqrt((pw2(x) + pw2(y))/pw2(r))*Ks23))/
            (sqrt(pw2(r))*(pw2(x) + pw2(y)));

          real_t K22 = (2*x*y*((pw2(x) + pw2(y))*Ks13 + 
                               sqrt((pw2(x) + pw2(y))/pw2(r))*z*Ks23))/
            (sqrt(pw2(r))*pow(pw2(x) + pw2(y),2));
          real_t K23 = (x*(z*Ks13 - sqrt((pw2(x) + pw2(y))/pw2(r))*Ks23))/
            (sqrt(pw2(r))*(pw2(x) + pw2(y)));
          real_t K33 = 0;

          
          
          real_t K = gammai11 * K11 + gammai22 * K22 + gammai33 * K33
            + 2.0 * (gammai12 * K12 + gammai13 * K13 + gammai23 * K23);

          DIFFK_a(i, j, k) = K;
          A11_a(i, j, k) = (K11 - DIFFgamma11_a(i, j, k) * K / 3.0) * pw2(chi_a(i, j, k));
          A12_a(i, j, k) = (K12 - DIFFgamma12_a(i, j, k) * K / 3.0) * pw2(chi_a(i, j, k));
          A13_a(i, j, k) = (K13 - DIFFgamma13_a(i, j, k) * K / 3.0) * pw2(chi_a(i, j, k));
          A22_a(i, j, k) = (K22 - DIFFgamma22_a(i, j, k) * K / 3.0) * pw2(chi_a(i, j, k));
          A23_a(i, j, k) = (K23 - DIFFgamma23_a(i, j, k) * K / 3.0) * pw2(chi_a(i, j, k));
          A33_a(i, j, k) = (K33 - DIFFgamma33_a(i, j, k) * K / 3.0) * pw2(chi_a(i, j, k));

          DIFFgamma11_a(i, j, k) *= pw2(chi_a(i, j, k));
          DIFFgamma22_a(i, j, k) *= pw2(chi_a(i, j, k));
          DIFFgamma33_a(i, j, k) *= pw2(chi_a(i, j, k));
          DIFFgamma12_a(i, j, k) *= pw2(chi_a(i, j, k));
          DIFFgamma13_a(i, j, k) *= pw2(chi_a(i, j, k));
          DIFFgamma23_a(i, j, k) *= pw2(chi_a(i, j, k));

          DIFFgamma11_a(i, j, k) -= 1.0;
          DIFFgamma22_a(i, j, k) -= 1.0;
          DIFFgamma33_a(i, j, k) -= 1.0;

          chi_a(i, j, k) -= 1.0;

        }
      }
    }
                
  }
  
}

void bssn_ic_static_blackhole(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
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

  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
       hierarchy->getGridGeometry()));
   TBOX_ASSERT(grid_geometry_);
   geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;


    const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

    
    std::shared_ptr<pdat::CellData<real_t> > chi_p_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFchi_p_idx)));

    const hier::Box& box = chi_p_pdata->getGhostBox();
    
    std::shared_ptr<pdat::CellData<real_t> > chi_a_pdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
         patch->getPatchData(DIFFchi_a_idx)));

    // std::shared_ptr<pdat::CellData<real_t> > chi_f_pdata(
    //   SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
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
            = pw2(2.0*norm / (1+ 2.0 * norm)) - 1.0;
        }
      }
    }
                
  }
  
}

void bssn_ic_ds_blackhole(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn)
{
# if ! USE_BSSN_SHIFT
  std::cerr << "Waning! USE_BSSN_SHIFT is suggested to be enabled in blackhole test " << std::endl;
# endif

  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  std::random_device rd;
  std::mt19937 gen(7.0 /*rd()*/);
  
  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
       hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;


     
  real_t H = 0.18, t0 = 0;
  real_t a = exp(H * t0), m = 1;


  idx_t DIFFchi_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFchi"), variable_db->getContext("ACTIVE"));
  idx_t DIFFchi_p_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFchi"), variable_db->getContext("PREVIOUS"));

  idx_t DIFFalpha_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFalpha"), variable_db->getContext("ACTIVE"));
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
    const std::shared_ptr<hier::Patch> & patch = *pit;

    std::shared_ptr<pdat::CellData<real_t> > DIFFchi_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFchi_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFchi_p_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFchi_p_idx)));

    std::shared_ptr<pdat::CellData<real_t> > DIFFalpha_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFalpha_a_idx)));

    const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry()));

    
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma12_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma12_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma13_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma13_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma22_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma23_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma23_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFgamma33_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFgamma33_a_idx)));

    std::shared_ptr<pdat::CellData<real_t> > A11_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A11_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A12_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A12_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A13_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A13_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A22_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A22_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A23_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A23_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > A33_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(A33_a_idx)));

    std::shared_ptr<pdat::CellData<real_t> > DIFFK_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFK_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > Gamma1_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(Gamma1_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > Gamma2_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(Gamma2_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > Gamma3_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(Gamma3_a_idx)));

    arr_t DIFFchi_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFchi_a_pdata->getArrayData());
    arr_t DIFFchi_p = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFchi_p_pdata->getArrayData());
        
    arr_t DIFFalpha_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFalpha_a_pdata->getArrayData());     

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


    const double * domain_lower = &grid_geometry.getXLower()[0];
    const double * domain_upper = &grid_geometry.getXUpper()[0];

    //const double *dx = &grid_geometry.getDx()[0];
    const double *dx = &patch_geom->getDx()[0];
      
    double L[3];


    for(int i = 0 ; i < 3; i++)
      L[i] = domain_upper[i] - domain_lower[i];
    
    
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
          real_t x = (dx[0] * ((real_t)i + 0.5)) - L[0] / 2.0 ;
          real_t y = (dx[1] * ((real_t)j + 0.5)) - L[1] / 2.0 ;
          real_t z = (dx[2] * ((real_t)k + 0.5)) - L[2] / 2.0 ;

          real_t r = sqrt(x*x + y*y + z*z);

          real_t ksi = m / (2.0 * a * r);

          //DIFFchi_p(i,j,k) = DIFFchi_a(i,j,k)
          //= pw2(2.0*r / (1+ 2.0 * r)) - 1.0;

          DIFFchi_p(i, j, k) = DIFFchi_a(i,j,k) = pw2((2.0*a*r)/(2.0*a*r + m)) / a - 1.0;
          DIFFK_a(i, j, k) = -3 * H;
          //          DIFFalpha_a(i, j, k) = - pw2( (1-ksi)/(1+ksi) ) - 1.0;
          bssn->K0 = -3 * H;
        }
      }
    }
                
  }
  
  
}

  
} // namespace cosmo
