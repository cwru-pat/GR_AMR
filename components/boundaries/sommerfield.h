
#include "../../cosmo_includes.h"
#include "../../cosmo_ps.h"

using namespace SAMRAI;

namespace cosmo{

class SommerfieldBD:CosmoPatchStrategy
{
 public:
  virtual void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double fill_time,
      const hier::IntVector& ghost_width_to_fill);
   hier::IntVector
   getRefineOpStencilWidth(
      const tbox::Dimension& dim) const;
   /* virtual void */
   /* preprocessRefineBoxes( */
   /*    hier::Patch& fine, */
   /*    const hier::Patch& coarse, */
   /*    const hier::BoxContainer& fine_boxes, */
   /*    const hier::IntVector& ratio); */
   virtual void
   preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);
   /* virtual void */
   /* postprocessRefineBoxes( */
   /*    hier::Patch& fine, */
   /*    const hier::Patch& coarse, */
   /*    const hier::BoxContainer& fine_boxes, */
   /*    const hier::IntVector& ratio); */
   virtual void
   postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);

   
   SommerfieldBD(
     const tbox::Dimension& dim,
     std::string object_name = std::string());
   vector<int> target_id_list;
   bool is_time_dependent;
};

}
