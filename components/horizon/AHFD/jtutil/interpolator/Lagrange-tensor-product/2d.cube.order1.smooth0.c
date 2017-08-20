/* $Header$ */

#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "../../util_ErrorCodes.h"
#include "../../../AHFD_macros.h"
#include "../InterpLocalUniform.h"
#include "../common/structs.h"
#include "../common/load.h"
#include "../common/evaluate.h"
#include "../common/store.h"

/* function prototype */
#define FUNCTION_NAME			AEILocalInterp_U_LagTP_2cube_10
#include "../template.h"

#define N_DIMS				2
#define MOLECULE_MIN_M			0
#define MOLECULE_MAX_M			1
#define MOLECULE_SIZE			2

/* which derivative ops do we support? */
#define HAVE_OP_I
#define HAVE_OP_DX
#define HAVE_OP_DY
/* n.b. no 2nd derivatives for linear interpolation! */

#define XYZ				x, y
#define FP_XYZ				fp x, fp y
#define STRIDE_IJK			stride_i, stride_j
#define JACOBIAN_MIJK_STRIDE		Jacobian_mi_stride, Jacobian_mj_stride

#define DATA_STRUCT			data_struct_2d_cube_size2
#define COEFFS_STRUCT			coeffs_struct_2d_cube_size2

#define LOAD_DATA_REAL			AEILocalInterp_load_2dcube2_r
#define LOAD_DATA_REAL4			AEILocalInterp_load_2dcube2_r4
#define LOAD_DATA_REAL8			AEILocalInterp_load_2dcube2_r8
#define LOAD_DATA_REAL16		AEILocalInterp_load_2dcube2_r16
#define LOAD_DATA_COMPLEX		AEILocalInterp_load_2dcube2_c
#define LOAD_DATA_COMPLEX8		AEILocalInterp_load_2dcube2_c8
#define LOAD_DATA_COMPLEX16		AEILocalInterp_load_2dcube2_c16
#define LOAD_DATA_COMPLEX32		AEILocalInterp_load_2dcube2_c32

#define EVALUATE_MOLECULE		AEILocalInterp_eval_2dcube2

#define STORE_COEFFS			AEILocalInterp_store_2dcube2

/* note pathnames are all relative to "../template.c" */
#define COEFFS_I_COMPUTE_FILE_NAME	"Lagrange-tensor-product/2d.coeffs/2d.cube.order1.smooth0/coeffs-I.compute.c"
#define COEFFS_DX_COMPUTE_FILE_NAME	"Lagrange-tensor-product/2d.coeffs/2d.cube.order1.smooth0/coeffs-dx.compute.c"
#define COEFFS_DY_COMPUTE_FILE_NAME	"Lagrange-tensor-product/2d.coeffs/2d.cube.order1.smooth0/coeffs-dy.compute.c"

/* actual code */
#include "../template.c"
