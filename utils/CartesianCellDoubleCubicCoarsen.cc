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
#include "CartesianCellDoubleCubicCoarsen.h"
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

CartesianCellDoubleCubicCoarsen::CartesianCellDoubleCubicCoarsen():
  hier::CoarsenOperator("CUBIC_COARSEN")
{

}


CartesianCellDoubleCubicCoarsen::~CartesianCellDoubleCubicCoarsen()
{
}

double CartesianCellDoubleCubicCoarsen::CINT(double u, double p0, double p1, double p2, double p3) const
{
  return 0.5*(
    (u*u*(2.0 - u) - u)*p0
    + (u*u*(3.0*u - 5.0) + 2)*p1
    + (u*u*(4.0 - 3.0*u) + u)*p2
    + u*u*(u - 1.0)*p3
  );
}

  
int
CartesianCellDoubleCubicCoarsen::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianCellDoubleCubicCoarsen::getStencilWidth(const tbox::Dimension& dim) const
{
  return hier::IntVector(dim, 1);
}

void
CartesianCellDoubleCubicCoarsen::coarsen(
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
           int il = i * ratio[0];
           int jl = j * ratio[1];
           int kl = k * ratio[2];

           double id = 0.5 ;
           double jd = 0.5 ;
           double kd = 0.5 ;
           
           double F_i_j_kd[16] = {0}, F_i_jd_kd[4] = {0};

           for(int si=0; si<4; ++si)
             for(int sj=0; sj<4; ++sj)
             {
               F_i_j_kd[si*4+sj] =
                 CINT(kd,
                      farray(il+si-1, jl+sj-1, kl-1), farray(il+si-1, jl+sj-1, kl+0),
                      farray(il+si-1, jl+sj-1, kl+1), farray(il+si-1, jl+sj-1, kl+2));
             }

           for(int si=0; si<4; ++si)
             F_i_jd_kd[si] = CINT(jd, F_i_j_kd[si*4+0], F_i_j_kd[si*4+1], F_i_j_kd[si*4+2], F_i_j_kd[si*4+3]);

           carray(i, j, k) = CINT(id, F_i_jd_kd[0], F_i_jd_kd[1], F_i_jd_kd[2], F_i_jd_kd[3]);

         }
       }
     }
   }
   else
   {
     TBOX_ERROR("CartesianCellDoubleCubicCoarsen error...\n"
                << "dim > 3 not supported." << std::endl);
   }


}

}
}
