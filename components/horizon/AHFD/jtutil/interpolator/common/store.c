/* store.c -- store molecule coefficients in Jacobian */
/* $Header$ */

#include "../../../AHFD_macros.h"
#include "../InterpLocalUniform.h"

#include "structs.h"
#include "store.h"

/******************************************************************************/

/*
 * 1-D routines
 */

#undef  COEFF
#define COEFF(mi)	Jacobian_ptr[Jacobian_mi_stride*mi]

void AEILocalInterp_store_1dcube2
      (fp factor, const struct coeffs_struct_1d_cube_size2 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride)
{
#include "1d.cube.size2/store-coeffs.c"
}

void AEILocalInterp_store_1dcube3
      (fp factor, const struct coeffs_struct_1d_cube_size3 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride)
{
#include "1d.cube.size3/store-coeffs.c"
}

void AEILocalInterp_store_1dcube4
      (fp factor, const struct coeffs_struct_1d_cube_size4 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride)
{
#include "1d.cube.size4/store-coeffs.c"
}

void AEILocalInterp_store_1dcube5
      (fp factor, const struct coeffs_struct_1d_cube_size5 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride)
{
#include "1d.cube.size5/store-coeffs.c"
}

void AEILocalInterp_store_1dcube6
      (fp factor, const struct coeffs_struct_1d_cube_size6 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride)
{
#include "1d.cube.size6/store-coeffs.c"
}

void AEILocalInterp_store_1dcube7
      (fp factor, const struct coeffs_struct_1d_cube_size7 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride)
{
#include "1d.cube.size7/store-coeffs.c"
}

/******************************************************************************/

/*
 * 2-D routines
 */

#undef  COEFF
#define COEFF(mi,mj)	Jacobian_ptr[   Jacobian_mi_stride*mi \
				      + Jacobian_mj_stride*mj ]

void AEILocalInterp_store_2dcube2
      (fp factor, const struct coeffs_struct_2d_cube_size2 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride)
{
#include "2d.cube.size2/store-coeffs.c"
}

void AEILocalInterp_store_2dcube3
      (fp factor, const struct coeffs_struct_2d_cube_size3 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride)
{
#include "2d.cube.size3/store-coeffs.c"
}

void AEILocalInterp_store_2dcube4
      (fp factor, const struct coeffs_struct_2d_cube_size4 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride)
{
#include "2d.cube.size4/store-coeffs.c"
}

void AEILocalInterp_store_2dcube5
      (fp factor, const struct coeffs_struct_2d_cube_size5 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride)
{
#include "2d.cube.size5/store-coeffs.c"
}

void AEILocalInterp_store_2dcube6
      (fp factor, const struct coeffs_struct_2d_cube_size6 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride)
{
#include "2d.cube.size6/store-coeffs.c"
}

/******************************************************************************/

/*
 * 3-D routines
 */

#undef  COEFF
#define COEFF(mi,mj,mk)	Jacobian_ptr[   Jacobian_mi_stride*mi \
				      + Jacobian_mj_stride*mj \
				      + Jacobian_mk_stride*mk ]

void AEILocalInterp_store_3dcube2
      (fp factor, const struct coeffs_struct_3d_cube_size2 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride, int Jacobian_mk_stride)
{
#include "3d.cube.size2/store-coeffs.c"
}

void AEILocalInterp_store_3dcube3
      (fp factor, const struct coeffs_struct_3d_cube_size3 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride, int Jacobian_mk_stride)
{
#include "3d.cube.size3/store-coeffs.c"
}

void AEILocalInterp_store_3dcube4
      (fp factor, const struct coeffs_struct_3d_cube_size4 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride, int Jacobian_mk_stride)
{
#include "3d.cube.size4/store-coeffs.c"
}

void AEILocalInterp_store_3dcube5
      (fp factor, const struct coeffs_struct_3d_cube_size5 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride, int Jacobian_mk_stride)
{
#include "3d.cube.size5/store-coeffs.c"
}

void AEILocalInterp_store_3dcube6
      (fp factor, const struct coeffs_struct_3d_cube_size6 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride, int Jacobian_mk_stride)
{
#include "3d.cube.size6/store-coeffs.c"
}
