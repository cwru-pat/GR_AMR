#include "bssn_ic.h"
#include "../../cosmo_includes.h"

using namespace SAMRAI;

namespace cosmo
{


void bssn_ic_static_blackhole(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
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

          chi_p(i,j,k) =chi_a(i,j,k)
            = 1/pw2((1.0 + 1.0/(2.0*norm))) - 1.0;
        }
      }
    }
                
  }
  
}
  
} // namespace cosmo
