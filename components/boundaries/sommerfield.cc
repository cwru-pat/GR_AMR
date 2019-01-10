#include "../../cosmo_includes.h"
#include "sommerfield.h"
#include "../bssn/bssn.h"

using namespace SAMRAI;

namespace cosmo{

/**
 * @brief   some comments for Sommerfield boundary:
 *          1. because it would be complicated to deal with inter patch communication 
 *             if we set the physical boundary of sommerfield on ghost cells, we use 
 *             most outer part of the inner cells of each patch on coarsest level as physical 
 *             boundaries. As the following example with ghost_width = 2, 
 *             where "*" represents cells on physical boundary, "x" represents
 *             the "real" inner cells. There are also ghost cells enclose all the "*" cells, but they are irrelevent 
 *             in our realization of Sommerfield boundary condition so are not shown here.
 *
 *            ***********
 *            ***********
 *            **xxxxxxx**
 *            **xxxxxxx**
 *            **xxxxxxx** 
 *            **xxxxxxx**
 *            **xxxxxxx**
 *            ***********
 *            ***********
 *
 *          2. this also means it would be unpredictable if we refine the most outer "*" cells,
 *             since the interpolation will be unprecise. So, let's try to avoid this case.
*/
SommerfieldBD::SommerfieldBD(
  const tbox::Dimension& dim_in,
  std::string object_name_in):
  CosmoPatchStrategy(dim_in, object_name_in)
   // if boundary is time dependent, it would not do anything when setting physical boundary conditions
{
  is_time_dependent = true;
}

SommerfieldBD::~SommerfieldBD() {
}

  
void SommerfieldBD::setPhysicalBoundaryConditions(
  hier::Patch& patch,
  const double fill_time,
  const hier::IntVector& ghost_width_to_fill)
{
  // do not do anything when setting physical boundaries
  // because they are evolved seperatedly in every RK step
  return;
}

/*
 * @brief should be useless
 */
void SommerfieldBD::preprocessRefine(
  hier::Patch& fine,
  const hier::Patch& coarse,
  const hier::Box& fine_box,
  const hier::IntVector& ratio)
{
  std::shared_ptr<hier::PatchGeometry> fine_geom (fine.getPatchGeometry());
  idx_t codim = 1;
  const std::vector<hier::BoundaryBox> & codim1_boxes =
    fine_geom->getCodimensionBoundaries(codim);
  const idx_t n_codim1_boxes = static_cast<idx_t>(codim1_boxes.size());

  //the return value of getBox() is const box, which can not be refined
  hier::Box coarse_box(coarse.getBox());

  coarse_box.refine(ratio);
  

  for(int bn = 0 ; bn < n_codim1_boxes; bn++)
  {
    const hier::Box & boundary_fill_box(
      fine_geom->getBoundaryFillBox(
        codim1_boxes[bn], fine_box * coarse_box, hier::IntVector(dim,STENCIL_ORDER)));

    if(boundary_fill_box.empty()) continue;
    
    const idx_t * lower = &boundary_fill_box.lower()[0];
    const idx_t * upper = &boundary_fill_box.upper()[0];

    //idx_t l_idx = boundary_fill_box.getLocationIndex();
    


    //initialize dx for each patch
    //const real_t * dx = &(fine_geom->getDx())[0];

    for (std::vector<int>::iterator it = target_id_list.begin() ;
         it != target_id_list.end(); ++it)
    {
      std::shared_ptr<pdat::CellData<double>> fine_pdata(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          fine.getPatchData(*it)));

      std::shared_ptr<pdat::CellData<double>> coarse_pdata(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          coarse.getPatchData(*it)));

      arr_t fine_a(pdat::ArrayDataAccess::access<DIM, double>(
                     fine_pdata->getArrayData()));

      arr_t coarse_a(pdat::ArrayDataAccess::access<DIM, double>(
                     coarse_pdata->getArrayData()));

      
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            fine_a(i,j,k) = coarse_a(i/2,j/2,k/2);
          }
        }
      }
    }
  }
  /******************codim = 2********************************/

  codim = 2;
  const std::vector<hier::BoundaryBox> & codim2_boxes =
    fine_geom->getCodimensionBoundaries(codim);
  const idx_t n_codim2_boxes = static_cast<idx_t>(codim2_boxes.size());


  for(int bn = 0 ; bn < n_codim2_boxes; bn++)
  {
    const hier::Box & boundary_fill_box(
      fine_geom->getBoundaryFillBox(
        codim2_boxes[bn], fine_box * coarse_box, hier::IntVector(dim,STENCIL_ORDER)));
    if(boundary_fill_box.empty()) continue;
    
    const idx_t * lower = &boundary_fill_box.lower()[0];
    const idx_t * upper = &boundary_fill_box.upper()[0];

    //    idx_t l_idx = codim2_boxes[bn].getLocationIndex();
    


    //initialize dx for each patch
    //const real_t * dx = &(fine_geom->getDx())[0];

    for (std::vector<int>::iterator it = target_id_list.begin() ;
         it != target_id_list.end(); ++it)
    {
      std::shared_ptr<pdat::CellData<double>> fine_pdata(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          fine.getPatchData(*it)));

      std::shared_ptr<pdat::CellData<double>> coarse_pdata(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          coarse.getPatchData(*it)));

      arr_t fine_a(pdat::ArrayDataAccess::access<DIM, double>(
                     fine_pdata->getArrayData()));

      arr_t coarse_a(pdat::ArrayDataAccess::access<DIM, double>(
                     coarse_pdata->getArrayData()));

      
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            fine_a(i,j,k) = coarse_a(i/2,j/2,k/2);
          }
        }
      }
    }
  }
  /********************************codim = 3************************/
  codim = 3;
  const std::vector<hier::BoundaryBox> & codim3_boxes =
    fine_geom->getCodimensionBoundaries(codim);
  const idx_t n_codim3_boxes = static_cast<idx_t>(codim3_boxes.size());


  for(int bn = 0 ; bn < n_codim3_boxes; bn++)
  {
    const hier::Box & boundary_fill_box( 
      fine_geom->getBoundaryFillBox(
        codim3_boxes[bn], fine_box * coarse_box, hier::IntVector(dim,STENCIL_ORDER)));

    if(boundary_fill_box.empty()) continue;
    
    const idx_t * lower = &boundary_fill_box.lower()[0];
    const idx_t * upper = &boundary_fill_box.upper()[0];

    //    idx_t l_idx = boundary_fill_box.getLocationIndex();
    


    //initialize dx for each patch
    //    const real_t * dx = &(fine_geom->getDx())[0];

    for (std::vector<int>::iterator it = target_id_list.begin() ;
         it != target_id_list.end(); ++it)
    {
      std::shared_ptr<pdat::CellData<double>> fine_pdata(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          fine.getPatchData(*it)));

      std::shared_ptr<pdat::CellData<double>> coarse_pdata(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          coarse.getPatchData(*it)));

      arr_t fine_a(pdat::ArrayDataAccess::access<DIM, double>(
                     fine_pdata->getArrayData()));

      arr_t coarse_a(pdat::ArrayDataAccess::access<DIM, double>(
                     coarse_pdata->getArrayData()));

      
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            fine_a(i,j,k) = coarse_a(i/2,j/2,k/2);
          }
        }
      }
    }
  }

}
void SommerfieldBD::postprocessRefine(
  hier::Patch& fine,
  const hier::Patch& coarse,
  const hier::Box& fine_box,
  const hier::IntVector& ratio)
{
  
}


hier::IntVector SommerfieldBD::getRefineOpStencilWidth(
  const tbox::Dimension& dim) const
{
  return hier::IntVector::getZero(dim);
}


}
