/* load.h -- load molecule-sized piece of input array into struct data */
/* $Header$ */

/*
 * prerequisite headers:
 *	"cctk.h"
 *	"../InterpLocalUniform.h"
 *	"structs.h"
 */

/******************************************************************************/

/*
 * 1-D load routines
 */

#undef  INT_STRIDE_IJK
#define INT_STRIDE_IJK			int stride_i

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_1dcube2_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_1d_cube_size2
#include "load-template.h"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_1dcube3_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_1d_cube_size3
#include "load-template.h"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_1dcube4_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_1d_cube_size4
#include "load-template.h"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_1dcube5_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_1d_cube_size5
#include "load-template.h"
#undef  LOAD_FUNCTION_NAME_PREFIX

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_1dcube6_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_1d_cube_size6
#include "load-template.h"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_1dcube7_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_1d_cube_size7
#include "load-template.h"

/******************************************************************************/

/*
 * 2-D load routines
 */

#undef  INT_STRIDE_IJK
#define INT_STRIDE_IJK			int stride_i, int stride_j

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_2dcube2_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_2d_cube_size2
#include "load-template.h"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_2dcube3_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_2d_cube_size3
#include "load-template.h"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_2dcube4_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_2d_cube_size4
#include "load-template.h"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_2dcube5_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_2d_cube_size5
#include "load-template.h"
#undef  LOAD_FUNCTION_NAME_PREFIX

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_2dcube6_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_2d_cube_size6
#include "load-template.h"

/******************************************************************************/

/*
 * 3-D load routines
 */

#undef  INT_STRIDE_IJK
#define INT_STRIDE_IJK			int stride_i, int stride_j, int stride_k

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_3dcube2_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_3d_cube_size2
#include "load-template.h"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_3dcube3_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_3d_cube_size3
#include "load-template.h"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_3dcube3_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_3d_cube_size3
#include "load-template.h"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_3dcube4_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_3d_cube_size4
#include "load-template.h"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_3dcube5_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_3d_cube_size5
#include "load-template.h"
#undef  LOAD_FUNCTION_NAME_PREFIX

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_3dcube6_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_3d_cube_size6
#include "load-template.h"

/******************************************************************************/

/*
 * We don't want to leave DATA_STRUCT defined -- this would confuse
 * later code that wants to define it and include "load-template.c"
 */
#undef  DATA_STRUCT
