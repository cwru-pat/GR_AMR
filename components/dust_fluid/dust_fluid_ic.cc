#include "dust_fluid.h"
#include "../../cosmo_includes.h"
#include "../../utils/math.h"


using namespace SAMRAI;

namespace cosmo
{

void DustFluid::dust_fluid_ic_set_fluid_for_BHL(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    idx_t ln, 
    std::shared_ptr<tbox::Database> cosmo_dust_fluid_db)
{
  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  double M = cosmo_dust_fluid_db->getDoubleWithDefault("M", 1);
  double a = cosmo_dust_fluid_db->getDoubleWithDefault("spin", 0.5);
  double K_c = cosmo_dust_fluid_db->getDoubleWithDefault("K_c", 1);
  double relaxation_tolerance = cosmo_dust_fluid_db->getDoubleWithDefault("relaxation_tolerance", 1e-8);
  int num_vcycles = cosmo_dust_fluid_db->getIntegerWithDefault("num_vcycles", 4);
  int max_depth = cosmo_dust_fluid_db->getIntegerWithDefault("max_depth", 4);

  double l = cosmo_dust_fluid_db->getDoubleWithDefault("l", 1);
  double sigma = cosmo_dust_fluid_db->getDoubleWithDefault("sigma", 3.5);

  if(!USE_BSSN_SHIFT)
    TBOX_ERROR("Must enable shift for blackhole simulation!\n");
  
  bssn_ic_kerr_BHL_CTT(hierarchy,ln, M, a, K_c, relaxation_tolerance
                       , num_vcycles, max_depth, l, sigma);

  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  idx_t DIFFchi_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFchi"), variable_db->getContext("ACTIVE"));;

  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    initPData(patch);
    initMDA(patch);
    
    const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

    std::shared_ptr<pdat::CellData<real_t> > chi_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFchi_a_idx)));
    arr_t DIFFchi_a =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        chi_a_pdata->getArrayData());

    const hier::Box& box = chi_a_pdata->getGhostBox();
    
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
          DF_D_a(i, j, k) = 1.0 ;
          DF_E_a(i, j, k) = 1.0 ;
          DF_S1_a(i, j, k) = 0;
          DF_S2_a(i, j, k) = 0;
          DF_S3_a(i, j, k) = 0;
          
        }
      }
    }
  }
  
}


}
