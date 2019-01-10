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
#include "CartesianCellDoubleCubicRefine.h"

#include <float.h>
#include <cmath>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/pdat/ArrayDataAccess.h"
#include "math.h"

/*
 *************************************************************************
 *
 * External declarations for FORTRAN  routines.
 *
 *************************************************************************
 */
extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

// in cartrefine3d.f:
void SAMRAI_F77_FUNC(cartquarefcelldoub3d, CARTQUAREFCELLDOUB3D) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
}

namespace SAMRAI {
namespace geom {

// using namespace std;
// tricubic interpolation
CartesianCellDoubleCubicRefine::CartesianCellDoubleCubicRefine():
  hier::RefineOperator("CUBIC_REFINE")
{

}

CartesianCellDoubleCubicRefine::~CartesianCellDoubleCubicRefine()
{
}

double CartesianCellDoubleCubicRefine::CINT(double u, double p0, double p1, double p2, double p3) const
{
  return 0.5*(
    (u*u*(2.0 - u) - u)*p0
    + (u*u*(3.0*u - 5.0) + 2)*p1
    + (u*u*(4.0 - 3.0*u) + u)*p2
    + u*u*(u - 1.0)*p3
  );
}
  
int
CartesianCellDoubleCubicRefine::getOperatorPriority() const
{
   return 0;
}

  
hier::IntVector
CartesianCellDoubleCubicRefine::getStencilWidth(const tbox::Dimension& dim) const
{
  return hier::IntVector(dim, 2);
}

void
CartesianCellDoubleCubicRefine::refine(
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
CartesianCellDoubleCubicRefine::refine(
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
   const int * ifirstc = &coarse_box.lower()[0];
   const int * ilastc = &coarse_box.upper()[0];
   const int * ifirstf = &fine_box.lower()[0];
   const int * ilastf = &fine_box.upper()[0];
   
   // initialize interpolation coefficient,
   // corresponds to u_{i-1}, u_{i} and u{i+1} respectively
   // only apply to the case
   // that interpolation ratio  = 2

   
   const double *cdx = &cgeom->getDx()[0];
   
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
           int il = (i < 0) ? ((i+1) / ratio[0] - 1): (i / ratio[0]);
           int jl = (j < 0) ? ((j+1) / ratio[1] - 1): (j / ratio[1]);
           int kl = (k < 0) ? ((k+1) / ratio[2] - 1): (k / ratio[2]);

           // get relative position to coarser points
           // will be 0 or 1 for left or right
           int pi  = 1 - i + il * ratio[0];
           int pj  = 1 - j + jl * ratio[1];
           int pk  = 1 - k + kl * ratio[2];

           il -= pi;
           jl -= pj;
           kl -= pk;
           
           double id = (pi == 1) ? (3.0 / 4.0) : (1.0 / 4.0);
           double jd = (pj == 1) ? (3.0 / 4.0) : (1.0 / 4.0);
           double kd = (pk == 1) ? (3.0 / 4.0) : (1.0 / 4.0);
                                   
           if( (pi != 0 && pi != 1) || (pj != 0 && pj != 1) || (pk != 0 && pk != 1))
             TBOX_ERROR("relative position is not 0 or 1 when doing quadratic interpolation");

           farray(i,j,k) = 0;
           double F_i_j_kd[16] = {0}, F_i_jd_kd[4] = {0};

           for(int si=0; si<4; ++si)
             for(int sj=0; sj<4; ++sj)
             {
               F_i_j_kd[si*4+sj] =
                 CINT(kd,
                      carray(il+si-1, jl+sj-1, kl-1), carray(il+si-1, jl+sj-1, kl+0),
                      carray(il+si-1, jl+sj-1, kl+1), carray(il+si-1, jl+sj-1, kl+2));
             }

           for(int si=0; si<4; ++si)
             F_i_jd_kd[si] = CINT(jd, F_i_j_kd[si*4+0], F_i_j_kd[si*4+1], F_i_j_kd[si*4+2], F_i_j_kd[si*4+3]);

           
           farray(i, j, k) = CINT(id, F_i_jd_kd[0], F_i_jd_kd[1], F_i_jd_kd[2], F_i_jd_kd[3]);
           
         }
       }
     }
   }
   else
   {
     TBOX_ERROR("CartesianCellDoubleCubicRefine error... \n"
                << "dim > 3 not supported." << std::endl);
   }

}

}
}
