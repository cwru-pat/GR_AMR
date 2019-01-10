#include "../../cosmo_includes.h"
#include "static_ic.h"

using namespace SAMRAI;

namespace cosmo
{

void dust_ic_set_random(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, std::shared_ptr<tbox::Database> cosmo_static_db)
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


  const double * domain_lower = &grid_geometry.getXLower()[0];
  const double * domain_upper = &grid_geometry.getXUpper()[0];

  real_t L[3];
  
  for(int i = 0 ; i < DIM; i++)
    L[i] = domain_upper[i] - domain_lower[i];

  
  ICsData icd = cosmo_get_ICsData(cosmo_static_db, L[0]);
  
  tbox::pout<<"Generating ICs with peak at k = "<<icd.peak_k<<"\n";
  tbox::pout<<"Generating ICs with peak amp. = "<<icd.peak_amplitude<<"\n";

  idx_t DIFFchi_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFchi"), variable_db->getContext("ACTIVE"));
  idx_t DIFFr_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFr"), variable_db->getContext("ACTIVE"));

  idx_t DIFFD_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFD"), variable_db->getContext("ACTIVE"));
  idx_t DIFFK_a_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("DIFFK"), variable_db->getContext("ACTIVE"));

  set_gaussian_random_field(
    hierarchy, ln, DIFFchi_a_idx, &icd);  


  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    std::shared_ptr<pdat::CellData<real_t> > DIFFchi_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFchi_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFr_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFr_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFD_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFD_a_idx)));
    std::shared_ptr<pdat::CellData<real_t> > DIFFK_a_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(DIFFK_a_idx)));


    arr_t DIFFchi_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFchi_a_pdata->getArrayData());     
    arr_t DIFFr_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFr_a_pdata->getArrayData());     
    arr_t DIFFD_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFD_a_pdata->getArrayData());     
    arr_t DIFFK_a = pdat::ArrayDataAccess::access<DIM, double>(
      DIFFK_a_pdata->getArrayData());     

  
    const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

    
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
          DIFFr_a(i,j,k)= -0.5/PI/(
            pow(1.0 + DIFFchi_a(i,j,k), 5.0)
          )*(
            double_derivative(i, j, k, 1, 1, DIFFchi_a,dx)
            + double_derivative(i, j, k, 2, 2, DIFFchi_a,dx)
            + double_derivative(i, j, k, 3, 3, DIFFchi_a,dx)
          );
          if(tbox::MathUtilities< double >::isNaN(DIFFr_a(i,j,k)))
          {
            TBOX_ERROR("Error: NaN DIFFr.");
          }

        }
      }
    }
    
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          DIFFchi_a(i,j,k) = 1.0 / pw2(1.0+DIFFchi_a(i,j,k)) - 1.0;
        }
      }
    }
    
    real_t min = icd.rho_K_matter;
    real_t max = min;
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          real_t rho_FRW = icd.rho_K_matter;
          real_t DIFFr = DIFFr_a(i,j,k);
          real_t rho = rho_FRW + DIFFr;
          // phi_FRW = 0
          real_t DIFFchi = DIFFchi_a(i,j,k);
          // phi = DIFFphi
          // DIFFK = 0

          DIFFD_a(i,j,k) =
            rho_FRW*(1.0 / pw3(DIFFchi + 1.0) - 1.0) + DIFFr/pw3(DIFFchi + 1.0);

          if(rho < min)
          {
            min = rho;
          }
          if(rho > max)
          {
            max = rho;
          }
          if(tbox::MathUtilities< double >::isNaN(rho))
          {
            TBOX_ERROR("Error: NaN energy density.");
          }
        }
      }
    }
    tbox::pout<< "Minimum fluid density: " <<min<<"\n";
    tbox::pout<< "Maximum fluid density: " <<max<<"\n";
    
    if(min < 0.0)
    {
      TBOX_ERROR("Error: negative density in some regions with min "<<min<<"\n");
    }
    
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {

          real_t rho_FRW = icd.rho_K_matter;
          real_t D_FRW = rho_FRW; // on initial slice

          DIFFr_a(i,j,k) += rho_FRW;

          DIFFK_a(i,j,k) = -sqrt(24.0*PI*rho_FRW);

          DIFFD_a(i,j,k) += D_FRW;
        }
      }
    }
  }
  
}
  
  
} // namespace cosmo
