/* InterpLocalUniform.h -- private stuff for this interpolator */
/* $Header$ */

/******************************************************************************/

/*
 * misc stuff that Jonathan Thornburg likes to use
 */

/* FIXME: C99 defines <stdbool.h>; we should really check if this has */
/*        already been included, and if so, not duplicate it */
//typedef int bool;
#include <stdbool.h>
/* make if-else symmetrical:  if (blah) then { ... } else { ... } */
#define then	/* empty */

/******************************************************************************/

/* number of integers in the range [x,y] inclusive */
#define HOW_MANY_IN_RANGE(x,y)  ((y) - (x) + 1)

/* is the integer x even/odd? */
#define IS_EVEN(x)	(((x) & 0x1) == 0)
#define IS_ODD (x)	(((x) & 0x1) != 0)

/* round floating-point value to nearest integer */
/* ... result is expressed as floating point! */
/* ... needs <math.h> for floor() */
#define ROUND_TO_INTEGER__RESULT_IS_FP(x)	floor((x) + 0.5)

/******************************************************************************/

/*
 * compile-time upper bounds for sizing arrays etc
 */

/*
 * We must have 1 <= N_dims <= MAX_N_DIMS.  Alas, there are lots of places
 * in "template.c" where code explicitly enumerates all possible N_DIMS
 * values at compile-time, so changing (eg raising) this limit requires
 * checking all preprocessor uses of N_DIMS as well as MAX_N_DIMS. :( :(
 */
#define MAX_N_DIMS	3

/* a "boundary" is the combination of a dimension and a min/max "side" */
#define MAX_N_BOUNDARIES	(2*MAX_N_DIMS)

/*
 * if molecule_family_string is a C string specifying a molecule family
 * (i.e. value associated with the "molecule_family" key in the parameter
 * table), we must have
 *	strlen(molecule_family_string) <= MAX_MOLECULE_FAMILY_STRLEN
 * n.b. exceeding this won't cause a buffer overflow, it will "just"
 *      cause the string to be truncated (and probably not recognized
 *      by the interpolator)
 */
#define MAX_MOLECULE_FAMILY_STRLEN	20

/*
 * N_MOLECULE_FAMILIES is the number of distinct molecule families
 */
enum	molecule_family
	{
	molecule_family_cube = 0,
	N_MOLECULE_FAMILIES	/* this must be the last entry in the enum */
	};

/*
 * We must have 1 <= order <= MAX_ORDER.
 */
#define MAX_ORDER	6

/*
 * We must have 0 <= smoothing <= MAX_SMOOTHING
 */
#define MAX_SMOOTHING	0

/******************************************************************************/

/*
 * other compile-time settings
 */

/* defaults for boundary_{off_centering,extrapolation}_tolerance[] */
#ifdef NOT_YET
  #define LAGRANGE_BNDRY_OFF_CNTR_TOL_DEF	999.0
  #define LAGRANGE_BNDRY_EXTRAP_TOL_DEF		1.0e-10
  #define HERMITE_BNDRY_OFF_CNTR_TOL_DEF	1.0e-10
  #define HERMITE_BNDRY_EXTRAP_TOL_DEF		0.0
#else
  #define LAGRANGE_BNDRY_OFF_CNTR_TOL_DEF	999.0
  #define LAGRANGE_BNDRY_EXTRAP_TOL_DEF		1.0e-10
  #define HERMITE_BNDRY_OFF_CNTR_TOL_DEF	999.0
  #define HERMITE_BNDRY_EXTRAP_TOL_DEF		1.0e-10
#endif

/* CCTK_VWarn() severity level for error/warning messages */
#define BUG_MSG_SEVERITY_LEVEL			0
#define ERROR_MSG_SEVERITY_LEVEL		0
#define WARNING_MSG_SEVERITY_LEVEL		1

/* CCTK_Abort() exit code for internal error (interpolator bug) aborts */
#define BUG_ABORT_CODE				42

/******************************************************************************/

/*
 * data structures used in multiple files
 */

/* number of real "parts" in a complex number */
#define COMPLEX_N_PARTS 2

#ifdef __cplusplus
extern "C" 
{
#endif

struct	error_info
	{
	/* did we find  error_point_status  in the parameter table? */
	bool found_per_point_status;

	/* should we print a Cactus level WARNING_MSG_SEVERITY_LEVEL */
	/* warning message if we find a point in error? */
	bool print_warning_msg;

	/* NULL pointer to skip per-point status, or */
	/* --> array of size N_interp_points to be set to per-point status */
	CCTK_INT* per_point_status;

	/* count of the number of points in error */
	CCTK_INT error_count;
	};

struct	molecule_structure_flags
	{
	bool MSS_is_fn_of_interp_coords;
	bool MSS_is_fn_of_which_operation;
	bool MSS_is_fn_of_input_array_values;
	bool Jacobian_is_fn_of_input_array_values;
	};

struct	molecule_min_max_m_info
	{
	CCTK_INT molecule_min_m[MAX_N_DIMS];
	CCTK_INT molecule_max_m[MAX_N_DIMS];
	};

struct	Jacobian_info
	{
	CCTK_REAL** Jacobian_pointer;
	CCTK_INT*   Jacobian_offset;
	CCTK_INT Jacobian_interp_point_stride;
	CCTK_INT Jacobian_m_strides[MAX_N_DIMS];
	CCTK_INT Jacobian_part_stride;
	};

/**************************************/

/*
 * error codes for AEILocalInterp_molecule_posn()
 */

/* x < minimum allowable x in grid */
#define MOLECULE_POSN_ERROR_X_LT_MIN	(-1)

/* x > maximum allowable x in grid */
#define MOLECULE_POSN_ERROR_X_GT_MAX	(-2)

/* grid is smaller than molecule */
#define MOLECULE_POSN_ERROR_GRID_TINY	(-3)

/* someone passed us a NaN or other non-finite floating-point number */
#define MOLECULE_POSN_ERROR_NAN		(-4)

/* grid spacing $\Delta x$ is zero in at least one axis */
#define MOLECULE_POSN_ERROR_DX_ZERO	(-5)

/******************************************************************************/

/*
 * prototypes for functions visible from multiple files
 */

/* functions in "startup.c" */
int AEILocalInterp_U_Startup(void);

/* functions in "InterpLocalUniform.c" */
/* ... these functions are registered by code in "startup.c", */
/*     then later called from the flesh */
int AEILocalInterp_U_Lagrange_TP(int N_dims,
				 int param_table_handle,
				 /***** coordinate system *****/
				 const CCTK_REAL coord_origin[],
				 const CCTK_REAL coord_delta[],
				 /***** interpolation points *****/
				 int N_interp_points,
				 int interp_coords_type_code,
				 const void *const interp_coords[],
				 /***** input arrays *****/
				 int N_input_arrays,
				 const CCTK_INT input_array_dims[],
				 const CCTK_INT input_array_type_codes[],
				 const void *const input_arrays[],
				 /***** output arrays *****/
				 int N_output_arrays,
				 const CCTK_INT output_array_type_codes[],
				 void *const output_arrays[]);
int AEILocalInterp_U_Lagrange_MD(int N_dims,
				 int param_table_handle,
				 /***** coordinate system *****/
				 const CCTK_REAL coord_origin[],
				 const CCTK_REAL coord_delta[],
				 /***** interpolation points *****/
				 int N_interp_points,
				 int interp_coords_type_code,
				 const void *const interp_coords[],
				 /***** input arrays *****/
				 int N_input_arrays,
				 const CCTK_INT input_array_dims[],
				 const CCTK_INT input_array_type_codes[],
				 const void *const input_arrays[],
				 /***** output arrays *****/
				 int N_output_arrays,
				 const CCTK_INT output_array_type_codes[],
				 void *const output_arrays[]);
int AEILocalInterp_U_Hermite(int N_dims,
			     int param_table_handle,
			     /***** coordinate system *****/
			     const CCTK_REAL coord_origin[],
			     const CCTK_REAL coord_delta[],
			     /***** interpolation points *****/
			     int N_interp_points,
			     int interp_coords_type_code,
			     const void *const interp_coords[],
			     /***** input arrays *****/
			     int N_input_arrays,
			     const CCTK_INT input_array_dims[],
			     const CCTK_INT input_array_type_codes[],
			     const void *const input_arrays[],
			     /***** output arrays *****/
			     int N_output_arrays,
			     const CCTK_INT output_array_type_codes[],
			     void *const output_arrays[]);

/* functions in "molecule_posn.c" */
int AEILocalInterp_molecule_posn(fp grid_origin, fp grid_delta,
				 int grid_i_min, int grid_i_max,
				 int molecule_size,
				 fp boundary_off_centering_tolerance_min,
				 fp boundary_off_centering_tolerance_max,
				 fp boundary_extrapolation_tolerance_min,
				 fp boundary_extrapolation_tolerance_max,
				 fp x,
				 int debug,
				 int* i_center, fp* x_rel);

/* functions in "util.c" */
int AEILocalInterp_decode_N_parts(int type_code);
int AEILocalInterp_get_int_param(const char* thorn_or_implementation_name,
				 const char* parameter_name);

#ifdef __cplusplus
}
#endif
