#ifndef COSMO_PS_H
#define COSMO_PS_H

#include "cosmo_includes.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

/*
 * Headers for basic SAMRAI objects used in this code.
 */

using namespace SAMRAI;

namespace cosmo
{
  
class CosmoPatchStrategy:public xfer::RefinePatchStrategy
{
 public:
   /*!
    * @brief Constructor using.
    *
    * Requires a concrete implementation of RobinBcCoefStrategy.
    *
    * @param dim
    * @param object_name Name of the object, for general referencing.
    * @param coef_strategy Coefficients strategy being helped.
    */

  
  explicit CosmoPatchStrategy(
    const tbox::Dimension& dim_in,
    std::string object_name_in = std::string());

  virtual ~CosmoPatchStrategy(
      void);
  

   void addTarget(idx_t idx);
   //@{ @name xfer::RefinePatchStrategy virtuals

   /* virtual void */
   /* setPhysicalBoundaryConditions( */
   /*    hier::Patch& patch, */
   /*    const double fill_time, */
   /*    const hier::IntVector& ghost_width_to_fill) = 0; */
   virtual hier::IntVector
   getRefineOpStencilWidth(
      const tbox::Dimension& dim) const = 0;
   /* virtual void */
   /* preprocessRefineBoxes( */
   /*    hier::Patch& fine, */
   /*    const hier::Patch& coarse, */
   /*    const hier::BoxContainer& fine_boxes, */
   /*    const hier::IntVector& ratio); */
   /* virtual void */
   /* preprocessRefine( */
   /*    hier::Patch& fine, */
   /*    const hier::Patch& coarse, */
   /*    const hier::Box& fine_box, */
   /*    const hier::IntVector& ratio); */
   /* virtual void */
   /* postprocessRefineBoxes( */
   /*    hier::Patch& fine, */
   /*    const hier::Patch& coarse, */
   /*    const hier::BoxContainer& fine_boxes, */
   /*    const hier::IntVector& ratio); */
   /* virtual void */
   /* postprocessRefine( */
   /*    hier::Patch& fine, */
   /*    const hier::Patch& coarse, */
   /*    const hier::Box& fine_box, */
   /*    const hier::IntVector& ratio); */

   std::vector<idx_t> target_id_list;


   const tbox::Dimension dim;
   std::string object_name;
   bool is_time_dependent;

};





}
#endif
