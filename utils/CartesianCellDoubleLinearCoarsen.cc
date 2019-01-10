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
#include "CartesianCellDoubleLinearCoarsen.h"
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

CartesianCellDoubleLinearCoarsen::CartesianCellDoubleLinearCoarsen():
  hier::CoarsenOperator("LINEAR_COARSEN")
{

}


CartesianCellDoubleLinearCoarsen::~CartesianCellDoubleLinearCoarsen()
{
}

int
CartesianCellDoubleLinearCoarsen::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianCellDoubleLinearCoarsen::getStencilWidth(const tbox::Dimension& dim) const
{
  return hier::IntVector(dim, 2);
}

void
CartesianCellDoubleLinearCoarsen::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(fine.getDim());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, coarse_box, ratio);

   std::shared_ptr<pdat::CellData<double> > cdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(src_component)));
   std::shared_ptr<pdat::CellData<double> > fdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(dst_component)));
   
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());

   MDA_Access<double, 3, MDA_OrderColMajor<3>> carray(
     pdat::ArrayDataAccess::access<3, double>(        
       cdata->getArrayData()));

   MDA_Access<double, 3, MDA_OrderColMajor<3>> farray(
     pdat::ArrayDataAccess::access<3, double>(        
       fdata->getArrayData()));

   
   const hier::Box fgbox(fdata->getGhostBox());

   const std::shared_ptr<CartesianPatchGeometry> cgeom(
      SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
         coarse.getPatchGeometry()));
   const std::shared_ptr<CartesianPatchGeometry> fgeom(
      SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
         fine.getPatchGeometry()));

   TBOX_ASSERT(cgeom);
   TBOX_ASSERT(fgeom);

   const hier::Box fine_box = hier::Box::refine(coarse_box, ratio);
   const int * ifirstc = &coarse_box.lower()[0];
   const int * ilastc = &coarse_box.upper()[0];
   const int * ifirstf = &fine_box.lower()[0];
   const int * ilastf = &fine_box.upper()[0];
   
   // initialize interpolation coefficient,
   // corresponds to u_{i-1}, u_{i} and u{i+1} respectively
   // only apply to the case
   // that interpolation ratio  = 2

   
   const double *fdx = &fgeom->getDx()[0];

   if ((dim == tbox::Dimension(3)))
   {
     #pragma omp parallel for collapse(2)
     for(int k = ifirstc[2]; k <= ilastc[2]; k++)
     {
       for(int j = ifirstc[1]; j <= ilastc[1]; j++)
       {
         for(int i = ifirstc[0]; i <= ilastc[0]; i++)
         {
           // generating corresponding coordinate on coarse mesh
           int fi = i * ratio[0];
           int fj = j * ratio[1];
           int fk = k * ratio[2];

           carray(i,j,k) = 0;

           // enumerating shift around coarse coordinate
           for(int sk = 0; sk <= 1; sk++)
           {
             double tempj = 0;
             for(int sj = 0; sj <= 1; sj++)
             {
               double tempi = 0;
               for(int si = 0; si <= 1; si++)
               {
                 tempi += farray(si + fi, sj + fj, sk + fk) * 0.5;
               }
               tempj += tempi * 0.5;
             }
             carray(i, j, k) += tempj * 0.5;
           }

         }
       }
     }
   }
   else
   {
     TBOX_ERROR("CartesianCellDoubleLinearRefine error...\n"
                << "dim > 3 not supported." << std::endl);
   }


}

}
}
