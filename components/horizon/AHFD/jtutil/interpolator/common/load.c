/* load.c -- load molecule-sized piece of input array into struct data */
/* $Header$ */

#include "../../../AHFD_macros.h"
#include "../InterpLocalUniform.h"

#include "structs.h"
#include "load.h"

/******************************************************************************/

/*
 * 1-D routines
 */

#undef  INT_STRIDE_IJK
#define INT_STRIDE_IJK		int stride_i
#undef  DATA_REAL
#define DATA_REAL(mi)		ptr[stride_i*mi]
#undef  DATA_COMPLEX
#define DATA_COMPLEX(mi)	ptr[stride_i*mi][part]

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_1dcube2_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_1d_cube_size2
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"1d.cube.size2/load-data.c"
#include "load-template.c"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_1dcube3_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_1d_cube_size3
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"1d.cube.size3/load-data.c"
#include "load-template.c"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_1dcube4_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_1d_cube_size4
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"1d.cube.size4/load-data.c"
#include "load-template.c"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_1dcube5_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_1d_cube_size5
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"1d.cube.size5/load-data.c"
#include "load-template.c"
#undef  LOAD_FUNCTION_NAME_PREFIX

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_1dcube6_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_1d_cube_size6
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"1d.cube.size6/load-data.c"
#include "load-template.c"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_1dcube7_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_1d_cube_size7
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"1d.cube.size7/load-data.c"
#include "load-template.c"

/******************************************************************************/

/*
 * 2-D routines
 */

#undef  INT_STRIDE_IJK
#define INT_STRIDE_IJK		int stride_i, int stride_j
#undef  DATA_REAL
#define DATA_REAL(mi,mj)	ptr[stride_i*mi + stride_j*mj]
#undef  DATA_COMPLEX
#define DATA_COMPLEX(mi,mj)	ptr[stride_i*mi + stride_j*mj][part]

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_2dcube2_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_2d_cube_size2
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"2d.cube.size2/load-data.c"
#include "load-template.c"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_2dcube3_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_2d_cube_size3
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"2d.cube.size3/load-data.c"
#include "load-template.c"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_2dcube4_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_2d_cube_size4
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"2d.cube.size4/load-data.c"
#include "load-template.c"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_2dcube5_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_2d_cube_size5
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"2d.cube.size5/load-data.c"
#include "load-template.c"
#undef  LOAD_FUNCTION_NAME_PREFIX

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_2dcube6_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_2d_cube_size6
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"2d.cube.size6/load-data.c"
#include "load-template.c"

/******************************************************************************/

/*
 * 3-D routines
 */

#undef  INT_STRIDE_IJK
#define INT_STRIDE_IJK		int stride_i, int stride_j, int stride_k
#undef  DATA_REAL
#define DATA_REAL(mi,mj,mk)	ptr[stride_i*mi + stride_j*mj + stride_k*mk]
#undef  DATA_COMPLEX
#define DATA_COMPLEX(mi,mj,mk)	ptr[stride_i*mi + stride_j*mj + stride_k*mk] \
				   [part]

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_3dcube2_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_3d_cube_size2
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"3d.cube.size2/load-data.c"
#include "load-template.c"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_3dcube3_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_3d_cube_size3
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"3d.cube.size3/load-data.c"
#include "load-template.c"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_3dcube4_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_3d_cube_size4
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"3d.cube.size4/load-data.c"
#include "load-template.c"

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_3dcube5_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_3d_cube_size5
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"3d.cube.size5/load-data.c"
#include "load-template.c"
#undef  LOAD_FUNCTION_NAME_PREFIX

#undef  LOAD_FUNCTION_NAME
#define LOAD_FUNCTION_NAME(type)	AEILocalInterp_load_3dcube6_ ## type
#undef  DATA_STRUCT
#define DATA_STRUCT			data_struct_3d_cube_size6
#undef  LOAD_DATA_FILE_NAME
#define LOAD_DATA_FILE_NAME		"3d.cube.size6/load-data.c"
#include "load-template.c"
