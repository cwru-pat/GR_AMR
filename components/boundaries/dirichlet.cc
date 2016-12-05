#include "cosmo_includes.h"
#include "../cosmo_bd.h"

/*
 * Headers for basic SAMRAI objects used in this code.
 */
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/TimeInterpolateOperator.h"
#include "SAMRAI/xfer/PatchLevelFillPattern.h"
#include "SAMRAI/xfer/PatchLevelBorderFillPattern.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"
#include "SAMRAI/xfer/PatchLevelInteriorFillPattern.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/geom/CartesianCellDoubleLinearRefine.h"

/*
 * Headers for major algorithm/data structure objects from SAMRAI
 */
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"

using namespace SAMRAI;

namespace cosmo{

void DirichletBD:setPhysicalBoundaryConditions(
  hier::Patch& patch,
  const double fill_time,
  const hier::IntVector& ghost_width_to_fill);
{
  boost::shared_ptr<pdat::CellData<double> > uval_data(
      patch.getPatchData(d_uval_data_id, getDataContext()),
      BOOST_CAST_TAG);
   const boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(                      
      patch.getPatchGeometry(),                               
      BOOST_CAST_TAG);                                         
   hier::IntVector ghost_cells(uval_data->getGhostCellWidth());

   // Set boundary conditions for cells corresponding to patch edges.
   const std::vector<hier::BoundaryBox>& edge_bdry =
     pgeom->getCodimensionBoundaries(Bdry::EDGE2D);
   
   for (int i = 0; i < static_cast<int>(edge_bdry.size()); i++) {
      // Boundary box       
      const hier::BoundaryBox &bbox = edge_bdry[i];
      // Utility for working with boundary boxes.
      hier::BoundaryBoxUtil bbox_util(bbox);
      // Compute box corresponding to bbox, but with correct width.
      hier::Box working_box(bbox.getDim());
      bbox_util.stretchBoxToGhostWidth(working_box,
                                       ghost_width_to_fill);
      // Fill data in the working box.
      uval_data->fillAll( 0.0, working_box );
   }
}

void DirichletBD:preprocessRefineBoxes(
  hier::Patch& fine,
  const hier::Patch& coarse,
  const hier::BoxContainer& fine_boxes,
  const hier::IntVector& ratio)
{
  
}
void DirichletBD:preprocessRefine(
  hier::Patch& fine,
  const hier::Patch& coarse,
  const hier::Box& fine_box,
  const hier::IntVector& ratio)
{
  
}
void DirichletBD:postprocessRefineBoxes(
  hier::Patch& fine,
  const hier::Patch& coarse,
  const hier::BoxContainer& fine_boxes,
  const hier::IntVector& ratio)
{
  
}
void DirichletBD:postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio)
{
  
}


}
