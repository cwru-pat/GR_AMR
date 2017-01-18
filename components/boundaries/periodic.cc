#include "../../cosmo_includes.h"
#include "periodic.h"
#include "../bssn/bssn.h"

using namespace SAMRAI;

namespace cosmo{

periodicBD::periodicBD(
  const tbox::Dimension& dim_in,
  std::string object_name_in):
  CosmoPatchStrategy(dim_in, object_name_in)
   // if boundary is time dependent, it would not do anything when setting physical boundary conditions
{
  is_time_dependent = false;
}

periodicBD::~periodicBD() {
}

  
void periodicBD::setPhysicalBoundaryConditions(
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
void periodicBD::preprocessRefine(
  hier::Patch& fine,
  const hier::Patch& coarse,
  const hier::Box& fine_box,
  const hier::IntVector& ratio)
{

}
void periodicBD::postprocessRefine(
  hier::Patch& fine,
  const hier::Patch& coarse,
  const hier::Box& fine_box,
  const hier::IntVector& ratio)
{
  
}


hier::IntVector periodicBD::getRefineOpStencilWidth(
  const tbox::Dimension& dim) const
{
  return hier::IntVector::getZero(dim);
}


}
