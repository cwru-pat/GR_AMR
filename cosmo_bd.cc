#ifndef COSMO_BD_H
#define COSMO_BD_H

#include "cosmo_includes.h"
#include "cosmo_bd.h"

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

namespace cosmo
{

CosmoBD: CosmoBD(
  const tbox::Dimension& dim_in):
  dim(dim_in)
{
  
}

CosmoBD: ~CosmoBD()
{
  
}

hier::IntVector CosmoBD:getRefineOpStencilWidth(
      const tbox::Dimension& dim) const;
{
  return hier::IntVector::getZero(dim);
}
  
}
