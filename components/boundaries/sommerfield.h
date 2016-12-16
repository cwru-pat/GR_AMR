#ifndef COSMO_SOMMERFIELD_H
#define COSMO_SOMMERFIELD_H

#include "../../cosmo_includes.h"
#include "../../cosmo_ps.h"

using namespace SAMRAI;

namespace cosmo{

class SommerfieldBD:public CosmoPatchStrategy
{
 public:

 SommerfieldBD(
    const tbox::Dimension& dim_in,
    std::string object_name);

  virtual ~SommerfieldBD(
    void);


  
  virtual void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double fill_time,
      const hier::IntVector& ghost_width_to_fill);

  virtual hier::IntVector
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

   
   std::vector<int> target_id_list;
   bool is_time_dependent;
};

}
#endif
