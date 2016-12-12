#ifndef COSMO_BD_H
#define COSMO_BD_H

#include "cosmo_includes.h"


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
  CosmoPatchStrategy();
   /*!
    * @brief Destructor.
    */
   virtual ~CosmoPatchStrategy(
      void);

   
   hier::IntVector getRefineOpStencilWidth(
      const tbox::Dimension& dim);

   void addTarget(idx_t idx);
   //@{ @name xfer::RefinePatchStrategy virtuals

   virtual void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double fill_time,
      const hier::IntVector& ghost_width_to_fill);
   hier::IntVector
   getRefineOpStencilWidth(
      const tbox::Dimension& dim) const;
   virtual void
   preprocessRefineBoxes(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::BoxContainer& fine_boxes,
      const hier::IntVector& ratio);
   virtual void
   preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);
   virtual void
   postprocessRefineBoxes(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::BoxContainer& fine_boxes,
      const hier::IntVector& ratio);
   virtual void
   postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);

   vector<idx_t> target_id_list;

   bool is_time_dependent;

   std::string object_name;
   const tbox::Dimension dim;
};





}
#endif
