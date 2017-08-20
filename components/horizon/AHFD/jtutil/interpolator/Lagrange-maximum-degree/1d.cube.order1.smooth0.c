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
#define FUNCTION_NAME			AEILocalInterp_U_LagMD_1cube_10
#include "../template.h"

#define N_DIMS				1
#define MOLECULE_MIN_M			0
#define MOLECULE_MAX_M			1
#define MOLECULE_SIZE			2

/* which derivative ops do we support? */
#define HAVE_OP_I
#define HAVE_OP_DX
/* n.b. no 2nd derivatives for linear interpolation! */

#define XYZ				x
#define FP_XYZ				fp x
#define STRIDE_IJK			stride_i
#define JACOBIAN_MIJK_STRIDE		Jacobian_mi_stride

#define DATA_STRUCT			data_struct_1d_cube_size2
#define COEFFS_STRUCT			coeffs_struct_1d_cube_size2

#define LOAD_DATA_REAL			AEILocalInterp_load_1dcube2_r
#define LOAD_DATA_REAL4			AEILocalInterp_load_1dcube2_r4
#define LOAD_DATA_REAL8			AEILocalInterp_load_1dcube2_r8
#define LOAD_DATA_REAL16		AEILocalInterp_load_1dcube2_r16
#define LOAD_DATA_COMPLEX		AEILocalInterp_load_1dcube2_c
#define LOAD_DATA_COMPLEX8		AEILocalInterp_load_1dcube2_c8
#define LOAD_DATA_COMPLEX16		AEILocalInterp_load_1dcube2_c16
#define LOAD_DATA_COMPLEX32		AEILocalInterp_load_1dcube2_c32

#define EVALUATE_MOLECULE		AEILocalInterp_eval_1dcube2

#define STORE_COEFFS			AEILocalInterp_store_1dcube2

/* note pathnames are all relative to "../template.c" */
#define COEFFS_I_COMPUTE_FILE_NAME	"Lagrange-maximum-degree/1d.coeffs/1d.cube.order1.smooth0/coeffs-I.compute.c"
#define COEFFS_DX_COMPUTE_FILE_NAME	"Lagrange-maximum-degree/1d.coeffs/1d.cube.order1.smooth0/coeffs-dx.compute.c"

/* actual code */
#include "../template.c"
