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
#include "CartesianCellDoubleQuadraticRefine.h"

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

CartesianCellDoubleQuadraticRefine::CartesianCellDoubleQuadraticRefine():
   hier::RefineOperator("QUADRATIC_REFINE")
{
}

CartesianCellDoubleQuadraticRefine::~CartesianCellDoubleQuadraticRefine()
{
}

int
CartesianCellDoubleQuadraticRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianCellDoubleQuadraticRefine::getStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getOne(dim);
}

void
CartesianCellDoubleQuadraticRefine::refine(
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
CartesianCellDoubleQuadraticRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(fine.getDim());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, fine_box, ratio);

   std::shared_ptr<pdat::CellData<double> > cdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(src_component)));
   std::shared_ptr<pdat::CellData<double> > fdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(dst_component)));
   
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());

   MDA_Access<double, 3, MDA_OrderColMajor<3>> carray(
     pdat::ArrayDataAccess::access<3, double>(        
       cdata->getArrayData()));

   MDA_Access<double, 3, MDA_OrderColMajor<3>> farray(
     pdat::ArrayDataAccess::access<3, double>(        
       fdata->getArrayData()));

   
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

   const int * ifirstf = &fine_box.lower()[0];
   const int * ilastf = &fine_box.upper()[0];

   // initialize interpolation coefficient,
   // corresponds to u_{i-1}, u_{i} and u{i+1} respectively
   // only apply to the case
   // that interpolation ratio  = 2

   
   const double coef[3][2] =
     {{5.0/32.0, -3.0/32.0},{15.0/16.0, 15.0/16.0},{-3.0/32.0, 5.0/32.0}};

   if ((dim == tbox::Dimension(3)))
   {
     #pragma omp parallel for collapse(2)
     for(int k = ifirstf[2]; k <= ilastf[2]; k++)
     {
       for(int j = ifirstf[1]; j <= ilastf[1]; j++)
       {
         for(int i = ifirstf[0]; i <= ilastf[0]; i++)
         {
           // generating corresponding coordinate on coarse mesh
           int ci = (i < 0) ? ((i+1) / ratio[0] - 1): (i / ratio[0]);
           int cj = (j < 0) ? ((j+1) / ratio[1] - 1): (j / ratio[1]);
           int ck = (k < 0) ? ((k+1) / ratio[2] - 1): (k / ratio[2]);

           // get relative position to coarser points
           // will be 0 or 1 for left or right
           int pi  = i - ci * ratio[0];
           int pj  = j - cj * ratio[1];
           int pk  = k - ck * ratio[2];

           if( (pi != 0 && pi != 1) || (pj != 0 && pj != 1) || (pk != 0 && pk != 1))
             TBOX_ERROR("relative position is not 0 or 1 when doing quadratic interpolation");

           farray(i,j,k) = 0;
           // enumerating shift around coarse coordinate
           for(int sk = -1; sk <= 1; sk++)
           {
             double tempj = 0;
             for(int sj = -1; sj <= 1; sj++)
             {
               double tempi = 0;
               for(int si = -1; si <= 1; si++)
               {
                 tempi += carray(ci+si, cj+sj,ck+sk) * coef[si + 1][pi];
               }
               tempj += tempi * coef[sj+1][pj];
             }
             farray(i,j,k) += tempj * coef[sk+1][pk];
           }
           
         }
       }
     }
   }
   else
   {
     TBOX_ERROR("CartesianCellDoubleQuadraticRefine error...\n"
                << "dim > 3 not supported." << std::endl);
   }
   //TBOX_ERROR("Terminating for test" << std::endl);

}

}
}
