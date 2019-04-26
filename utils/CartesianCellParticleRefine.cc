/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Linear refine operator for cell-centered double data on
 *                a Cartesian mesh.
 *
 ************************************************************************/
#include "CartesianCellParticleRefine.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/pdat/ArrayDataAccess.h"
/*
 *************************************************************************
 *
 * External declarations for FORTRAN  routines.
 *
 *************************************************************************
 */

namespace SAMRAI {
namespace geom {

// using namespace std;

CartesianCellParticleRefine::CartesianCellParticleRefine():
   hier::RefineOperator("PARTICLE_REFINE")
{
}

CartesianCellParticleRefine::~CartesianCellParticleRefine()
{
}

int
CartesianCellParticleRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianCellParticleRefine::getStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getOne(dim);
}

void
CartesianCellParticleRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::BoxOverlap& fine_overlap,
   const hier::IntVector& ratio) const
{
   const pdat::CellOverlap* t_overlap =
      CPP_CAST<const pdat::CellOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      refine(fine,
         coarse,
         dst_component,
         src_component,
         *b,
         ratio);
   }
}

void
CartesianCellParticleRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(fine.getDim());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, fine_box, ratio);

   std::shared_ptr<
     pdat::IndexData<ParticleContainer,
                     pdat::CellGeometry> > cdata(
                       SAMRAI_SHARED_PTR_CAST<pdat::IndexData<ParticleContainer, pdat::CellGeometry>,
                       hier::PatchData>(
                         coarse.getPatchData(src_component)));

   std::shared_ptr<
     pdat::IndexData<ParticleContainer,
                     pdat::CellGeometry> > fdata(
                       SAMRAI_SHARED_PTR_CAST<pdat::IndexData<ParticleContainer, pdat::CellGeometry>,
                       hier::PatchData>(
                         fine.getPatchData(dst_component)));
   
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(fdata);

   
   const hier::Box cgbox(cdata->getGhostBox());

   const std::shared_ptr<CartesianPatchGeometry> cgeom(
      SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
         coarse.getPatchGeometry()));
   const std::shared_ptr<CartesianPatchGeometry> fgeom(
      SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
         fine.getPatchGeometry()));

   TBOX_ASSERT(cgeom);
   TBOX_ASSERT(fgeom);

   const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);

   //   tbox::pout<<coarse_box<<" "<<fine_box<<"\n";
   pdat::CellIterator icend(pdat::CellGeometry::end(coarse_box));
   for (pdat::CellIterator ic(pdat::CellGeometry::begin(coarse_box));
                              ic != icend; ++ic)
   {
     // get pointer to this particle at Index
     ParticleContainer * p = cdata->getItem(*ic);

     // no particles here
     if( p == NULL ) continue;

     const double* cdx = cgeom->getDx();
     // calculating cell center coordinate
     double center_x = ( (double)((*ic)[0]) + 0.5) * cdx[0];
     double center_y = ( (double)((*ic)[1]) + 0.5) * cdx[1];
     double center_z = ( (double)((*ic)[2]) + 0.5) * cdx[2];

     // calculating corresponding idx on left bottom conner
     hier::Index f_idx(0, 0, 0);

     f_idx[0] = ((*ic)[0]*2 );
     f_idx[1] = ((*ic)[1]*2 );
     f_idx[2] = ((*ic)[2]*2 );

     
     // go through all particles
     for(std::list<RKParticle>::iterator it=p->p_list.begin();
         it != p->p_list.end(); ++it)
     {
       // calculating which finer cell the particle live in
       // by comparing its location of cell center
       int si = SIGN((*it).x_a[0] - center_x) > 0 ? 1 : 0;
       int sj = SIGN((*it).x_a[1] - center_y) > 0 ? 1 : 0;
       int sk = SIGN((*it).x_a[2] - center_z) > 0 ? 1 : 0;
       hier::Index shift(si,sj,sk);

       hier::Index p_idx = f_idx + shift;
       
       ParticleContainer *fp = fdata->getItem(p_idx);

       if(fp)
       {
         fp->addParticle(*it);
       }
       else
       {
         fp = new ParticleContainer(p_idx);
         fp->addParticle(*it);
         fdata->appendItem(p_idx, *fp);
       }
     }
     

   }

}

}
}
