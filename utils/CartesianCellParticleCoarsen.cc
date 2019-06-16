/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for cell-centered double data on
 *                a Cartesian mesh.
 *
 ************************************************************************/
#include "CartesianCellParticleCoarsen.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/pdat/ArrayDataAccess.h"

#include <float.h>
#include <cmath>
#include "math.h"
/*
 *************************************************************************
 *
 * External declarations for FORTRAN  routines.
 *
 *************************************************************************
 */


namespace SAMRAI {
namespace geom {

CartesianCellParticleCoarsen::CartesianCellParticleCoarsen():
  hier::CoarsenOperator("PARTICLE_COARSEN")
{
}


CartesianCellParticleCoarsen::~CartesianCellParticleCoarsen()
{
}

int
CartesianCellParticleCoarsen::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianCellParticleCoarsen::getStencilWidth(const tbox::Dimension& dim) const
{
  return hier::IntVector::getZero(dim);
}

void
CartesianCellParticleCoarsen::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(fine.getDim());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, coarse_box, ratio);

   std::shared_ptr<
     pdat::IndexData<ParticleContainer,
                     pdat::CellGeometry> > cdata(
                       SAMRAI_SHARED_PTR_CAST<pdat::IndexData<ParticleContainer, pdat::CellGeometry>,
                       hier::PatchData>(
                         coarse.getPatchData(dst_component)));

   std::shared_ptr<
     pdat::IndexData<ParticleContainer,
                     pdat::CellGeometry> > fdata(
                       SAMRAI_SHARED_PTR_CAST<pdat::IndexData<ParticleContainer, pdat::CellGeometry>,
                       hier::PatchData>(
                         fine.getPatchData(src_component)));
   
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata);

   
   const std::shared_ptr<CartesianPatchGeometry> cgeom(
      SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
         coarse.getPatchGeometry()));
   const std::shared_ptr<CartesianPatchGeometry> fgeom(
      SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
         fine.getPatchGeometry()));

   TBOX_ASSERT(cgeom);
   TBOX_ASSERT(fgeom);

   
   hier::Box fine_box = fdata->getGhostBox();
   hier::Box c_box = cdata->getGhostBox();

   // if(c_box.getGlobalId().getOwnerRank() == 1) return;
   
   pdat::CellIterator fiend(pdat::CellGeometry::end(fine_box));
   for (pdat::CellIterator fi(pdat::CellGeometry::begin(fine_box));
                              fi != fiend; ++fi)
   {
     // get pointer to this particle at Index
     ParticleContainer * p = fdata->getItem(*fi);

     if( p == NULL ) continue;

     // calculating corresponding idx on left bottom conner
     hier::Index c_idx(0, 0, 0);

     c_idx[0] = ((*fi)[0] >= 0) ? ((*fi)[0]/2 ) : ( ((*fi)[0] -1)/2);
     c_idx[1] = ((*fi)[1] >= 0) ? ((*fi)[1]/2 ) : ( ((*fi)[1] -1)/2);
     c_idx[2] = ((*fi)[2] >= 0) ? ((*fi)[2]/2 ) : ( ((*fi)[2] -1)/2);


     if(!c_box.contains(c_idx))
       TBOX_ERROR("Coarse index is not contained by patch data\n");

     // go through all particles
     for(std::list<RKParticle>::iterator it=p->p_list.begin();
         it != p->p_list.end(); ++it)
     {
       ParticleContainer *cp = cdata->getItem(c_idx);

       if(cp)
       {
         cp->addParticle(*it);
       }
       else
       {
         //std::cout<<"hhhhh"<<c_idx<<" "<<coarse.getBox()<<" "<<coarse_box<<"\n"<<std::flush;       
         cp = new ParticleContainer(c_idx);
         cp->addParticle(*it);
         cdata->appendItem(c_idx, *cp);
       }

     }
     
   }
   //   std::cout<<"Coarsen data gets "<<cdata->getNumberOfItems()<<" particles\n";
}

}
}
