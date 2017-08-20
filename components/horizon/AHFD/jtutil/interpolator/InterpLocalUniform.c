/* InterpLocalUniform.c -- driver for this interpolator */
/* $Header$ */
/*
** *****data structures and functions local to this file *****
**
 * AEILocalInterp_U_Lagrange_TP - driver for Lagrange interpolator (tensor prod)
 * AEILocalInterp_U_Lagrange_MD - driver for Lagrange interpolator (max degree)
 * AEILocalInterp_U_Hermite - driver for Hermite interpolator
**
** InterpLocalUniform - main driver routine
**
** check_boundary_tolerances - check boundary tolerances for consistency
** get_error_point_info - get per-point error reporting info from par table
** get_and_decode_molecule_family - get & decode molecule_family from par table
** get_molecule_positions - get molecule_positions from parameter table
** get_Jacobian_info - get Jacobian-query info from parameter table
** set_error_info - set error information in parameter table
** set_molecule_structure - set molecule structure info in parameter table
** set_molecule_min_max_m - set molecule size info in parameter table
**
** get_and_check_INT - get and range-check CCTK_INT from parameter table
** get_INT_array - get CCTK_INT array from parameter table
** get_REAL_array - get CCTK_REAL array from parameter table
 */

/*@@
  @file      InterpLocalUniform.c
  @date      22 Oct 2001
  @author    Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc
	Generalized Interpolation of processor-local arrays.

	This code interpolates a set of functions defined on a
	uniform N-dimensional grid, to a set of interpolation points.
	It can also do differentiation and/or smoothing in the process
	of the interpolation.

	Conceptually, the generalized interpolation is done by
	least-squares fitting a polynomial of the specified order
	to the data values, applying the differentiation (if any)
	to it, then evaluating the result at the interpolation points.
	Since this ultimately yields _some_ linear combination of the
	data values, we precompute the coefficients for efficiency.
	This is done with separate Maple-generated formulas for each
	combination of number of interpolation operator, number of
	dimensions, molecule family, order, amount of Savitzky-Golay
	smoothing to be done, and differentiation operator.
  @enddesc
  @version   $Id$
  @@*/

#include <stdio.h>
#include <stdlib.h>			/* malloc(), free() */
#include <string.h>
#include <stdbool.h>

#include "../../AHFD_macros.h"
#include "../util_String.h"
#include "../util_Table.h"
#include "InterpLocalUniform.h"

#include "Lagrange-tensor-product/all_prototypes.h"
#include "Lagrange-maximum-degree/all_prototypes.h"
#include "Hermite/all_prototypes.h"

/* the rcs ID and its dummy function to use it */
/* static const char* rcsid = "$Header$"; */
/* CCTK_FILEVERSION(AEITHorns_AEILocalInterp_src_InterpLocalUniform_c) */

/******************************************************************************/

/*
 * ***** data structures and functions local to this file *****
 */

/**************************************/

/*
 * data structures local to this file
 */

/* this enum describes which interpolation operator we're using */
enum	interp_operator
	{
	interp_operator_Lagrange_TP,	/* "Lagrange polynomial interpolation */
					/*  (tensor product)" */
	interp_operator_Lagrange_MD,	/* "Lagrange polynomial interpolation */
					/*  (maximum degree)" */
	interp_operator_Hermite,	/* "Hermite polynomial interpolation" */
	N_INTERP_OPERATORS	/* this must be the last entry in the enum */
	};

/**************************************/

/*
 * prototypes for functions local to this file
 */
static
  int InterpLocalUniform(enum interp_operator interp_operator,
			 const char* interp_operator_string,
			 CCTK_REAL boundary_off_cntr_tol_default,
			 CCTK_REAL boundary_extrap_tol_default,
		       /***** misc arguments *****/
			 int N_dims,
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

static
  void check_boundary_tolerances
	(int N_boundaries,
	 const CCTK_REAL boundary_off_centering_tolerance[MAX_N_BOUNDARIES],
	 const CCTK_REAL boundary_extrapolation_tolerance[MAX_N_BOUNDARIES]);
static
  int get_error_point_info(int param_table_handle,
			   struct error_info *p_error_info);
static
  int get_and_decode_molecule_family
	(int param_table_handle,
	 int buffer_size, char molecule_family_string_buffer[],
	 enum molecule_family *p_molecule_family);
static
  int get_molecule_positions
	(int param_table_handle,
	 int N_dims, CCTK_INT* molecule_positions_array[MAX_N_DIMS]);
static
  int get_Jacobian_info(int param_table_handle,
			int N_dims, int N_output_arrays,
			struct Jacobian_info* p_Jacobian_info);
static
  int set_error_info(int param_table_handle,
		     struct error_info* p_error_info);
static
  int set_molecule_structure
	(int param_table_handle,
	 const struct molecule_structure_flags* p_molecule_structure_flags);
static
  int set_molecule_min_max_m
	(int param_table_handle,
	 int N_dims,
	 const struct molecule_min_max_m_info* p_molecule_min_max_m_info);

static
  int get_and_check_INT(int param_table_handle, const char name[],
			bool mandatory_flag, int default_value,
			bool check_range_flag, int min_value, int max_value,
			const char max_value_string[],
			CCTK_INT* p_value);
static
  int get_INT_array(int param_table_handle, const char name[],
		    bool default_flag, int default_value,
		    int N, CCTK_INT buffer[],
		    bool* p_value_not_set);
static
  int get_REAL_array(int param_table_handle, const char name[],
		      CCTK_REAL default_value,
		      int N, CCTK_REAL buffer[]);

/**************************************/

/*
 * table of function pointers pointing to the actual interpolation functions
 */

/*
 * typedef  p_interp_fn_t  as a function pointer pointing to an
 * individual interpolator function of the sort defined by template.[ch]
 */
#undef  FUNCTION_NAME
#define FUNCTION_NAME	(*p_interp_fn_t)
typedef
  #include "template.h"

/* NULL (function) pointer of this type */
#define NULL_INTERP_FN_PTR        ((p_interp_fn_t) NULL)

/*
 * For each axis of this array which is marked with a "see above" comment
 * in the array declaration, the array size is declared as X+1, so the
 * legal subscripts are of course 0...X.  But for these axes we actually
 * only use 1...X.  (Having the array oversize like this slightly simplifies
 * the array subscripting.)
 *
 * In the comments on the initializers, "n.i." = "not implemented".
 */
static const p_interp_fn_t p_interp_fn_table[N_INTERP_OPERATORS]
					    [MAX_N_DIMS+1]	/* see above */
					    [N_MOLECULE_FAMILIES]
					    [MAX_ORDER+1]	/* see above */
					    [MAX_SMOOTHING+1]
    = {

	{
	/* interp_operator == interp_operator_Lagrange_TP */
	  {
	  /* N_dims = 0 ==> unused */
	    {
	    /* molecule_family = molecule_family_cube */
	      { NULL_INTERP_FN_PTR,                },	/* order=0 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=1 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=2 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=3 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=4 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=5 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=6 (unused) */
	    },
	  },

	  {
	  /* N_dims = 1 */
	    {
	    /* molecule_family = molecule_family_cube */
	      { NULL_INTERP_FN_PTR,                },	/* order=0 (unused) */
	      { AEILocalInterp_U_LagTP_1cube_10,   },	/* order=1 */
	      { AEILocalInterp_U_LagTP_1cube_20,   },	/* order=2 */
	      { AEILocalInterp_U_LagTP_1cube_30,   },	/* order=3 */
	      { AEILocalInterp_U_LagTP_1cube_40,   },	/* order=4 */
	      { AEILocalInterp_U_LagTP_1cube_50,   },	/* order=5 */
	      { AEILocalInterp_U_LagTP_1cube_60,   },	/* order=6 */
	    },
	  },

	  {
	  /* N_dims = 2 */
	    {
	    /* molecule_family = molecule_family_cube */
	      { NULL_INTERP_FN_PTR,                },	/* order=0 (unused) */
	      { AEILocalInterp_U_LagTP_2cube_10,   },	/* order=1 */
	      { AEILocalInterp_U_LagTP_2cube_20,   },	/* order=2 */
	      { AEILocalInterp_U_LagTP_2cube_30,   },	/* order=3 */
	      { AEILocalInterp_U_LagTP_2cube_40,   },	/* order=4 */
	      { NULL_INTERP_FN_PTR,                },	/* order=5 (n.i.) */
	      { NULL_INTERP_FN_PTR,                },	/* order=6 (n.i.) */
	    },
	  },

	  {
	  /* N_dims = 3 */
	    {
	    /* molecule_family = molecule_family_cube */
	      { NULL_INTERP_FN_PTR,                },	/* order=0 (unused) */
	      { AEILocalInterp_U_LagTP_3cube_10,   },	/* order=1 */
	      { AEILocalInterp_U_LagTP_3cube_20,   },	/* order=2 */
	      { AEILocalInterp_U_LagTP_3cube_30,   },	/* order=3 */
	      { AEILocalInterp_U_LagTP_3cube_40,   },	/* order=4 */
	      { NULL_INTERP_FN_PTR,                },	/* order=5 (n.i.) */
	      { NULL_INTERP_FN_PTR,                },	/* order=6 (n.i.) */
	    },
	  },

	/* end of interp_operator == interp_operator_Lagrange_TP */
	},

	{
	/* interp_operator == interp_operator_Lagrange_MD */
	  {
	  /* N_dims = 0 ==> unused */
	    {
	    /* molecule_family = molecule_family_cube */
	      { NULL_INTERP_FN_PTR,                },	/* order=0 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=1 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=2 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=3 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=4 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=5 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=6 (unused) */
	    },
	  },

	  {
	  /* N_dims = 1 */
	    {
	    /* molecule_family = molecule_family_cube */
	      { NULL_INTERP_FN_PTR,                },	/* order=0 (unused) */
	      { AEILocalInterp_U_LagMD_1cube_10,   },	/* order=1 */
	      { AEILocalInterp_U_LagMD_1cube_20,   },	/* order=2 */
	      { AEILocalInterp_U_LagMD_1cube_30,   },	/* order=3 */
	      { AEILocalInterp_U_LagMD_1cube_40,   },	/* order=4 */
	      { AEILocalInterp_U_LagMD_1cube_50,   },	/* order=5 */
	      { AEILocalInterp_U_LagMD_1cube_60,   },	/* order=6 */
	    },
	  },

	  {
	  /* N_dims = 2 */
	    {
	    /* molecule_family = molecule_family_cube */
	      { NULL_INTERP_FN_PTR,                },	/* order=0 (unused) */
	      { AEILocalInterp_U_LagMD_2cube_10,   },	/* order=1 */
	      { AEILocalInterp_U_LagMD_2cube_20,   },	/* order=2 */
	      { AEILocalInterp_U_LagMD_2cube_30,   },	/* order=3 */
	      { AEILocalInterp_U_LagMD_2cube_40,   },	/* order=4 */
	      { NULL_INTERP_FN_PTR,                },	/* order=5 (n.i.) */
	      { NULL_INTERP_FN_PTR,                },	/* order=6 (n.i.) */
	    },
	  },

	  {
	  /* N_dims = 3 */
	    {
	    /* molecule_family = molecule_family_cube */
	      { NULL_INTERP_FN_PTR,                },	/* order=0 (unused) */
	      { AEILocalInterp_U_LagMD_3cube_10,   },	/* order=1 */
	      { AEILocalInterp_U_LagMD_3cube_20,   },	/* order=2 */
	      { AEILocalInterp_U_LagMD_3cube_30,   },	/* order=3 */
	      { AEILocalInterp_U_LagMD_3cube_40,   },	/* order=4 */
	      { NULL_INTERP_FN_PTR,                },	/* order=5 (n.i.) */
	      { NULL_INTERP_FN_PTR,                },	/* order=6 (n.i.) */
	    },
	  },

	/* end of interp_operator == interp_operator_Lagrange_MD */
	},

	{
	/* interp_operator == interp_operator_Hermite */
	  {
	  /* N_dims = 0 ==> unused */
	    {
	    /* molecule_family = molecule_family_cube */
	      { NULL_INTERP_FN_PTR,                },	/* order=0 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=1 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=2 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=3 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=4 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=5 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=6 (unused) */
	    },
	  },

	  {
	  /* N_dims = 1 */
	    {
	    /* molecule_family = molecule_family_cube */
	      { NULL_INTERP_FN_PTR,                },	/* order=0 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=1 (unused) */
	      { AEILocalInterp_U_Herm_1cube_2,     },	/* order=2 */
	      { AEILocalInterp_U_Herm_1cube_3,     },	/* order=3 */
	      { AEILocalInterp_U_Herm_1cube_4,     },	/* order=4 */
	      { NULL_INTERP_FN_PTR,                },	/* order=5 (n.i.) */
	      { NULL_INTERP_FN_PTR,                },	/* order=6 (n.i.) */
	    },
	  },

	  {
	  /* N_dims = 2 */
	    {
	    /* molecule_family = molecule_family_cube */
	      { NULL_INTERP_FN_PTR,                },	/* order=0 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=1 (unused) */
	      { AEILocalInterp_U_Herm_2cube_2,     },	/* order=2 */
	      { AEILocalInterp_U_Herm_2cube_3,     },	/* order=3 */
	      { NULL_INTERP_FN_PTR,                },	/* order=4 (n.i.) */
	      { NULL_INTERP_FN_PTR,                },	/* order=5 (n.i.) */
	      { NULL_INTERP_FN_PTR,                },	/* order=6 (n.i.) */
	    },
	  },

	  {
	  /* N_dims = 3 */
	    {
	    /* molecule_family = molecule_family_cube */
	      { NULL_INTERP_FN_PTR,                },	/* order=0 (unused) */
	      { NULL_INTERP_FN_PTR,                },	/* order=1 (unused) */
	      { AEILocalInterp_U_Herm_3cube_2,     },	/* order=2 */
	      { AEILocalInterp_U_Herm_3cube_3,     },	/* order=3 */
	      { NULL_INTERP_FN_PTR,                },	/* order=4 (n.i.) */
	      { NULL_INTERP_FN_PTR,                },	/* order=5 (n.i.) */
	      { NULL_INTERP_FN_PTR,                },	/* order=6 (n.i.) */
	    },
	  },

	/* end of interp_operator == interp_operator_Hermite */
	},

      };

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*@@
  @routine    AEILocalInterp_U_Lagrange_TP
  @date       22 Oct 2001
  @author     Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc
	This function is a top-level driver for
	LocalInterp::CCTK_InterpLocalUniform().

	It does Lagrange interpolation of a set of input arrays defined
	on a uniform N-dimensional tensor-product grid, to a set of
	interpolation points.  For multidimensional interpolation, it
	uses the tensor-product basis for the interpolating function,
	so the interpolant is guaranteed to be continuous and to pass
	through the input data at the grid points (up to floating-point
	roundoff errors).

	This interpolator can also
	* do differentiation and/or smoothing in the process of the
	  interpolation
	* return information about the Jacobian of the interpolated
	  values with respect to the input arrays

	This function does nothing except pass all its arguments
	down to  InterpLocalUniform()  with an extra 2 arguments
	added to indicate the operator handle.  See that function
	for details on this function's arguments/results.
  @enddesc
  @@*/
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
				 void *const output_arrays[])
{
return InterpLocalUniform(interp_operator_Lagrange_TP,
			  "Lagrange polynomial interpolation (tensor product)",
			  LAGRANGE_BNDRY_OFF_CNTR_TOL_DEF,
			  LAGRANGE_BNDRY_EXTRAP_TOL_DEF,
			  N_dims,
			  param_table_handle,
			  coord_origin, coord_delta,
			  N_interp_points,
			    interp_coords_type_code,
			    interp_coords,
			  N_input_arrays,
			    input_array_dims,
			    input_array_type_codes,
			    input_arrays,
			  N_output_arrays,
			    output_array_type_codes,
			    output_arrays);
}

/******************************************************************************/

/*@@
  @routine    AEILocalInterp_U_Lagrange_MD
  @date       22 Oct 2001
  @author     Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc
	This function is a top-level driver for
	LocalInterp::CCTK_InterpLocalUniform().

	It does Lagrange interpolation of a set of input arrays defined
	on a uniform N-dimensional tensor-product grid, to a set of
	interpolation points.  For multidimensional interpolation, it
	uses the maximum-degree basis for the interpolating function.
	This means the interpolant is in general *not* continuous and
	in general does *not* pass through the input data at the grid
	points.

	This interpolator can also
	* do differentiation and/or smoothing in the process of the
	  interpolation
	* return information about the Jacobian of the interpolated
	  values with respect to the input arrays

	This function does nothing except pass all its arguments
	down to  InterpLocalUniform()  with an extra 2 arguments
	added to indicate the operator handle.  See that function
	for details on this function's arguments/results.
  @enddesc
  @@*/
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
				void *const output_arrays[])
{
return InterpLocalUniform(interp_operator_Lagrange_MD,
			  "Lagrange polynomial interpolation (maximum degree)",
			  LAGRANGE_BNDRY_OFF_CNTR_TOL_DEF,
			  LAGRANGE_BNDRY_EXTRAP_TOL_DEF,
			  N_dims,
			  param_table_handle,
			  coord_origin, coord_delta,
			  N_interp_points,
			    interp_coords_type_code,
			    interp_coords,
			  N_input_arrays,
			    input_array_dims,
			    input_array_type_codes,
			    input_arrays,
			  N_output_arrays,
			    output_array_type_codes,
			    output_arrays);
}

/******************************************************************************/

/*@@
  @routine    AEILocalInterp_U_Hermite
  @date       22 Oct 2001
  @author     Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc
	This function is a top-level driver for
	LocalInterp::CCTK_InterpLocalUniform().

	It does Hermite interpolation of a set of input arrays defined
	on a uniform N-dimensional tensor-product grid, to a set of
	interpolation points.  It can also do differentiation and/or
	smoothing in the process of the interpolation.  It can also
	return information about the Jacobian of the interpolated
	values with respect to the input arrays.

	This function does nothing except pass all its arguments
	down to  InterpLocalUniform()  with an extra 2 arguments
	added to indicate the operator handle.  See that function
	for details on this function's arguments/results.
  @enddesc
  @@*/
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
			    void *const output_arrays[])
{
return InterpLocalUniform(interp_operator_Hermite,
			  "Hermite polynomial interpolation",
			  HERMITE_BNDRY_OFF_CNTR_TOL_DEF,
			  HERMITE_BNDRY_EXTRAP_TOL_DEF,
			  N_dims,
			  param_table_handle,
			  coord_origin, coord_delta,
			  N_interp_points,
			    interp_coords_type_code,
			    interp_coords,
			  N_input_arrays,
			    input_array_dims,
			    input_array_type_codes,
			    input_arrays,
			  N_output_arrays,
			    output_array_type_codes,
			    output_arrays);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*@@
  @routine    InterpLocalUniform
  @date       22 Oct 2001
  @author     Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc
	This function is the main driver for
	LocalInterp::CCTK_InterpLocalUniform().

	It interpolates a set of input arrays defined on a uniform
	N-dimensional tensor-product grid, to a set of interpolation
	points.  It can also do differentiation and/or smoothing
	in the process of the interpolation.  It can also return
	information about the Jacobian of the interpolated values
	with respect to the input arrays.

	This function is just a driver: it validates arguments,
	extracts optional arguments from the parameter table,
	then decodes
	   (interp_operator, N_dims, molecule_family, order, smoothing)
	and calls the appropriate subfunction to do the actual
	interpolation, then finally stores some results back in
	the parameter table.
  @enddesc

  @hdate    28.Jan.2003
  @hauthor  Jonathan Thornburg <jthorn@aei.mpg.de>
  @hdesc    Support parameter-table entries
            @var N_boundary_points_to_omit,
            @var boundary_off_centering_tolerance, and
            @var boundary_extrapolation_tolerance.
  @endhdesc

  @hdate    1.Feb.2003
  @hauthor  Jonathan Thornburg <jthorn@aei.mpg.de>
  @hdesc    Change defaults to "no off-centering" for Hermite interpolator,
            split off most of this (huge) function into subfunctions
  @endhdesc

  ***** operator arguments *****

  @var		interp_operator
  @vdesc	describes the interpolation operator
  @vtype	enum interp_operator interp_operator
  @endvar

  @var		interp_operator_string
  @vdesc	gives the character-string name of the interpolation operator
		(this is used only for printing error messages)
  @vtype	const char* interp_operator_string
  @endvar

  @var		boundary_off_cntr_tol_default
  @vdesc	the default value for each element of
		@var boundary_off_centering_default[]
		if this key isn't present in the parameter table
  @vtype	CCTK_REAL boundary_off_cntr_tol_default
  @endvar

  @var		boundary_extrap_tol_default
  @vdesc	the default value for each element of
		@var boundary_extrapolation_default[]
		if this key isn't present in the parameter table
  @vtype	CCTK_REAL boundary_extrap_tol_default
  @endvar

  ***** misc arguments *****

  @var		N_dims
  @vdesc	dimensionality of the interpolation
  @vtype	int N_dims			(must be >= 1)
  @endvar

  @var		param_table_handle
  @vdesc	handle to a key-value table giving additonal parameters
		for the interpolation
  @vtype	int param_table_handle		(must be >= 0)
  @endvar

  ***** arguments specifying the coordinate system *****

  @var		coord_origin
  @vdesc	(pointer to) array[N_dims] of values giving the
		x,y,z,... coordinates which correspond to the
		integer input-array subscripts (0,0,0,...,0)
		(note there is no implication here that such a
		grid point need actually exist; the arrays just
		give the coordinates it would have if it did exist)
  @vtype	const CCTK_REAL coord_origin[N_dims]
  @endvar

  @var		coord_delta
  @vdesc	(pointer to) array[N_dims] of values giving the
		coordinate spacing of the grid
  @vtype	const CCTK_REAL coord_delta[N_dims]
  @endvar

  ***** arguments specifying the interpolation points *****

  @var		N_interp_points
  @vdesc number of interpolation points
  @vtype	int N_interp_points		(must be >= 0)
  @endvar

  @var		interp_coords_type_code
  @vdesc	one of the CCTK_VARIABLE_* codes giving the data
		type of the arrays pointed to by  interp_coords[]
  @vtype	int
  @endvar

  @var		interp_coords
  @vdesc	(pointer to) array[N_dims] of pointers
		to arrays[N_interp_points] giving
		x,y,z,... coordinates of interpolation points
  @vtype	const void *const interp_coords[N_dims]
  @endvar

  ***** arguments specifying the inputs (the data to be interpolated) *****

  @var		N_input_arrays
  @vdesc	number of arrays input to the interpolation
  @vtype	int N_input_arrays		(must be >= 0)
  @endvar

  @var		input_array_dims
  @vdesc	dimensions of the input arrays: unless overridden by
		entries in the parameter table, all input arrays are
		taken to have these dimensions, with [0] the most contiguous
		axis and [N_dims-1] the least contiguous axis, and
		array subscripts in the range 0 <= subscript < dims[axis]
  @vtype	const int input_array_dims[N_dims]
  @endvar

  @var		input_array_type_codes
  @vdesc	(pointer to) array of CCTK_VARIABLE_* codes
		giving the data types of the input arrays
  @vtype	const int input_array_type_codes[N_input_arrays]
  @endvar

  @var		input_arrays
  @vdesc	(pointer to) array[N_input_arrays] of pointers to input arrays
  @vtype	const void *const input_arrays[N_input_arrays]
  @endvar

  ***** arguments specifying the outputs (the interpolation results) *****

  @var		N_output_arrays
  @vdesc	number of arrays output from the interpolation
  @vtype	int N_output_arrays		(must be >= 0)
  @endvar

  @var		output_array_type_codes
  @vdesc	(pointer to) array of CCTK_VARIABLE_* codes
		giving the data types of the output arrays
  @vtype	const int output_array_type_codes[N_output_arrays]
  @endvar

  @var		output_arrays
  @vdesc	(pointer to) array[N_output_arrays] of pointers to output arrays
  @vtype	void *const output_arrays[N_output_arrays]
  @endvar

  ***** input arguments passed in the parameter table *****

  @var		order
  @vdesc	Sets the order of the interpolating polynomial
		(1=linear, 2=quadratic, 3=cubic, ...).
		This table entry is mandatory; all others are optional.
  @vtype	int order
  @endvar

  @var		N_boundary_points_to_omit
  @vdesc	If this key is present, then this array specifies how many
		input grid points to omit (not use for the interpolation)
		on each grid boundary.  The array elements are matched up
		with the axes and minimum/maximum ends of the grid in the
		order [x_min, x_max, y_min, y_max, z_min, z_max, ...].
		If this key isn't present, it defaults to all zeros, i.e.
		to use all the input grid points in the interpolation.
  @vtype	const CCTK_INT N_boundary_points_to_omit[2*N_dims]
  @endvar

  @var		boundary_off_centering_tolerance
  @vdesc	If this key is present, then the interpolator allows
		interpolation points to be up to (<=) this many grid
		spacings outside the default-centering region on each
		grid boundary, off-centering the interpolation molecules
		as necessary.
		The array elements are matched up with the axes and
		minimum/maximum ends of the grid in the order
		[x_min, x_max, y_min, y_max, z_min, z_max, ...].
		If this key isn't present, it defaults to all infinities,
		i.e. to place no restriction on the interpolation point.
  @vtype	const CCTK_REAL boundary_off_centering_tolerance[2*N_dims]
  @endvar

  @var		boundary_extrapolation_tolerance
  @vdesc	If this key is present, then the interpolator allows
		interpolation points to be up to (<=) this many grid
		spacings outside the valid-data region on each grid
		boundary, doing extrapolation instead of interpolation
		as necessary.
		The array elements are matched up with the axes and
		minimum/maximum ends of the grid in the order
		[x_min, x_max, y_min, y_max, z_min, z_max, ...].
		If this key isn't present, it defaults to all 1.0e-10,
		i.e. to allow up to 1e-10 grid spacings of extrapolation.
  @vtype	const CCTK_REAL boundary_extrapolation_tolerance[2*N_dims]
  @endvar

  @var		input_array_offsets
  @vdesc	If this key is present, the value should be a pointer to
		an array giving an "offset" for each input array, use in
		the subscripting computation:  for input array number in,
		this computation (given for 3D; the generalization to other
		numbers of dimensions is obvious) is
		  input_arrays[in][offset + i*istride + j*jstride + k*kstride]
		where
		  offset = input_array_offsets[in]
			   or is 0 if input_array_offsets is absent
		  (istride,jstride,kstride,...) = input_array_stride[]
		and where (i,j,k,...) run from input_array_min_subscripts[]
		to input_array_max_subscripts[] inclusive.
  @vtype	const CCTK_INT input_array_offsets[N_input_arrays]
  @endvar

  @var		input_array_strides
  @vdesc	(pointer to) array giving strides of input arrays
		(this is shared by all input arrays)
  @vtype	const CCTK_INT input_array_strides[N_dims]
  @endvar

  @var		input_array_min_subscripts
  @vdesc	(pointer to) array giving minimum subscripts of input arrays
		(this is shared by all input arrays)
  @vtype	const CCTK_INT input_array_min_subscripts[N_dims]
  @endvar

  @var		input_array_max_subscripts
  @vdesc	(pointer to) array giving maximum subscripts of input arrays
		(this is shared by all input arrays)
  @vtype	const CCTK_INT input_array_max_subscripts[N_dims]
  @endvar

  @var		operand_indices
  @vdesc	(pointer to) array of integer operand indices specifying
		which input array (0-origin indexing into  input_arrays[])
		is to be generalized-interpolated to produce this output
  @vtype	const CCTK_INT operand_indices[N_output_arrays]
  @endvar

  @var		operation_codes
  @vdesc	(pointer to) array of integer operation codes specifying
		what (if any) derivatives are to be takin in the interpolation
		process: 0=no derivative, 1=d/dx, 2=d/dy, 3=d/dz (for coords
		(x,y,z)), 11=d^2/dx^2, 12=21=d^2/dxdy, etc etc.
  @vtype	const CCTK_INT operation_codes[N_output_arrays]
  @endvar

  ***** output arguments passed in the parameter table *****

  @var		out_of_range_pt
  @vdesc	If an out-of-range point is detected, this table entry
		is set to its  pt  value.
  @vtype	CCTK_INT out_of_range_pt
  @vio		out
  @endvar

  @var		out_of_range_axis
  @vdesc	If an out-of-range point is detected, this table entry
		is set to the  axis  value which is out of range.
  @vtype	CCTK_INT out_of_range_axis
  @vio		out
  @endvar

  @var		out_of_range_end
  @vdesc	If an out-of-range point is detected, this table entry
		is set to -1/+1 if the point is out of range on the min/max
		end of the axis.
  @vtype	CCTK_INT out_of_range_end
  @vio		out
  @endvar

  @var		MSS_is_fn_of_interp_coords
  @vdesc	If the interpolation molecule size and/or shape vary
		with the interpolation coordinates, this table entry
		is set to 1.  Otherwise (i.e. if the interpolation
		molecule size and shape are independent of the
		interpolation coordinates), it's set to 0.
  @vtype	CCTK_INT MSS_is_fn_of_interp_coords
  @vio		out
  @endvar

  @var		MSS_is_fn_of_which_operation
  @vdesc	If the interpolator supports computing derivatives,
		and the interpolation molecule size and/or shape vary
		from one  operation_code[] value to another, this table
		entry is set to 1.  Otherwise (i.e. if the interpolator
		doesn't support computing derivatives, or if the
		interpolator does support computing derivatives but
		the interpolation molecule size and shape are independent   
		of the  operation_code[]  values), it's set to 0.
		Note that this query tests whether the molecule size
		and/or shape depend on  operation_codes[]  in general,
		independent of whether there are in fact any distinct values
		(or even any values at all) passed in  operation_codes[]
		in this particular interpolator call.  In other words,
		this is a query about the basic design of this
		interpolator, not about this particular call.
  @vtype	CCTK_INT MSS_is_fn_of_which_operation
  @vio		out
  @endvar

  @var		MSS_is_fn_of_input_array_values
  @vdesc	If the interpolation molecule size and/or shape vary
		with the actual floating-point values of the input
		arrays, this table entry is set to 1.  Otherwise
		(i.e. if the interpolation molecule size and shape
		are independent of the input array values; this is
		a necessary, but not sufficient, condition for the
		interpolation to be linear), it's set to 0.
  @vtype	CCTK_INT MSS_is_fn_of_input_array_values
  @vio		out
  @endvar

  @var		Jacobian_is_fn_of_input_array_values
  @vdesc	If the actual floating-point values of the Jacobian (*)
		(see the discussion of  Jacobian_pointer  below)
		depends on the actual floating-point values of the
		input arrays (i.e. if the interpolation is nonlinear),
		this table entry is set to 1.  Otherwise (i.e. if the
		interpolation is linear), it's set to 0.
  @vtype	CCTK_INT Jacobian_is_fn_of_input_array_values
  @vio		out
  @endvar

  @var		molecule_min_m
  @vdesc	If this key is present (the value doesn't matter),
  		then the value is (re)set to an array of N_dims integers
		giving the minimum molecule m coordinate in each axis
  @vtype	CCTK_INT molecule_min_m[N_dims]
  @vio		out
  @endvar

  @var		molecule_max_m
  @vdesc	If this key is present (the value doesn't matter),
  		then the value is (re)set to an array of N_dims integers
		giving the maximum molecule m coordinate in each axis
  @vtype	CCTK_INT molecule_max_m[N_dims]
  @vio		out
  @endvar

  @var		molecule_positions
  @vdesc	If this key is present, then the value should be an
		array of N_dims pointers to (caller-supplied) arrays
		of N_interp_points CCTK_INTs each; the interpolator
		will store the molecule positions in the pointed-to
		arrays.
  @vtype	CCTK_INT *const molecule_positions[N_dims]
  @vio		out
  @endvar

  @var		Jacobian_pointer
  @vdesc	If this key is present, then the value should be an
		array of N_output_arrays pointers to (caller-supplied)
		arrays of CCTK_INTs; the interpolator will store the
		Jacobian elements in the pointed-to arrays.  For output
		array number out, the subscripting computation (given
		for 3D; the generalization to other numbers of dimensions
		is obvious) is
		 Jacobian_pointer[out][offset
				       + pt*Jacobian_interp_point_stride
				       + mi*stride_i + mj*stride_j + mk*stride_k
				       + part*Jacobian_part_stride]
		where
		  offset = Jacobian_offset[out]
			   or is 0 if Jacobian_offset[] is absent
		  pt = the point number
		  (stride_i,stride_j,stride_k) = Jacobian_m_strides[]
		  part = 0 for real values,
			 0 for the real parts of complex values,
			 1 for the imaginary parts of complex values, or
			 0 if Jacobian_part_stride is absent
  @vtype	CCTK_REAL *const Jacobian_pointer[N_output_arrays]
  @endvar

  @var		Jacobian_offset
  @vdesc	If this key is present, then the value should be a
		pointer to an array of N_output_arrays CCTK_INTs giving
		offsets to use in the Jacobian subscripting computation.
		If this key is absent, the effect is as if a default
		array of all 0s had been supplied.
  @vtype	const CCTK_INT Jacobian_offset[N_output_arrays]
  @endvar

  @var		Jacobian_interp_point_stride
  @vdesc	The pt-stride for the Jacobian subscripting computation.
  @vtype	const CCTK_INT Jacobian_interp_point_stride
  @endvar

  @var		Jacobian_m_strides
  @vdesc	An array of N_dims CCTK_INTs giving the m strides along
		each axis for the Jacobian subscripting computation
  @vtype	const CCTK_INT Jacobian_m_strides[N_dims]
  @endvar

  @var		Jacobian_part_stride
  @vdesc	If this key is present, it gives the part stride for the
		Jacobian subscripting computation.  If this key is absent,
		the default value 0 is used (n.b. this is suitable only
		for real data arrays)
  @vtype	const CCTK_INT Jacobian_m_strides[N_dims]
  @endvar

  ***** return result *****

  @returntype	int
  @returndesc	0			successful interpolation
		UTIL_ERROR_BAD_INPUT	one of the input arguments is invalid
		UTIL_ERROR_NO_MEMORY	unable to malloc() temporary memory
		UTIL_ERROR_BAD_HANDLE	invalid parameter table handle
		CCTK_ERROR_INTERP_POINT_OUTSIDE
					interpolation point is out of range
					(in this case additional information
					is reported through the parameter table)
		or any error code returned by one of the Util_TableGet*()
		functions called by this function
  @endreturndesc
  @@*/
static
  int InterpLocalUniform(enum interp_operator interp_operator,
			 const char* interp_operator_string,
			 CCTK_REAL boundary_off_cntr_tol_default,
			 CCTK_REAL boundary_extrap_tol_default,
		       /***** misc arguments *****/
			 int N_dims,
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
			 void *const output_arrays[])
{
/*
 * Implementation Note:
 *
 * We malloc() several scratch arrays, some with sizes determined by
 * N_{input,output}_arrays.  Thus if N_{input,output}_arrays == 0, with
 * the obvious code we would malloc(0).  Alas, the C standard permits
 * malloc(0) to return a NULL pointer, which the usual malloc() idiom
 *	CCTK_INT *const p = malloc(N * sizeof(CCTK_INT));
 *	if (p == NULL)
 *	   then return UTIL_ERROR_NO_MEMORY
 * would falsely detect as an out-of-memory condition.
 *
 * As a workaround, we pad all our malloc request sizes, i.e.
 *	CCTK_INT* const p = (CCTK_INT*) malloc((N+1) * sizeof(CCTK_INT));
 *	if (p == NULL)
 *	   then return UTIL_ERROR_NO_MEMORY
 * This is a kludge, but so are the other plausible solutions. :( :(
 */
int N_input_arrays1  = N_input_arrays  + 1;
int N_output_arrays1 = N_output_arrays + 1;

int status, status1, status2;
bool value_not_set;

/******************************************************************************/

/*
 * basic sanity checks on input parameters
 */
if (    (N_dims <= 0)
     || (param_table_handle < 0)
     || ((N_interp_points > 0) && (coord_origin == NULL))
     || ((N_interp_points > 0) && (coord_delta == NULL))
     || (N_interp_points < 0)
     || ((N_interp_points > 0) && (interp_coords == NULL))
     || (N_input_arrays < 0)
	/* no check here on input_array_dims, */
	/* since it may be NULL if overridden by parameter-table stuff */
     || ((N_input_arrays > 0) && (input_array_type_codes == NULL))
     || ((N_input_arrays > 0) && (input_arrays == NULL))
     || (N_output_arrays < 0)
     || ((N_output_arrays > 0) && (output_array_type_codes == NULL))
     || ((N_output_arrays > 0) && (output_arrays == NULL))    )
   then {
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"CCTK_InterpLocalUniform(): invalid arguments\n"
"   (N_dims <= 0, param_table_handle < 0, N_input_arrays < 0,\n"
"    N_output_arrays < 0, and/or NULL pointers-that-shouldn't-be-NULL)!");
	return UTIL_ERROR_BAD_INPUT;			/*** ERROR RETURN ***/
	}

if (N_dims > MAX_N_DIMS)
   then {
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): implementation restriction: N_dims=%d\n"
"                              but we only support N_dims <= MAX_N_DIMS=%d!",
		   N_dims,
		   MAX_N_DIMS);
	return UTIL_ERROR_BAD_INPUT;			/*** ERROR RETURN ***/
	}
   
/******************************************************************************/

/*
 * if logging is requested,
 *    open a logging file if we haven't already done so
 */
  {
/******************************************************************************/

/*
 * get the parameters from the parameter table, filling in defaults as needed
 */

CCTK_INT debug;
status = get_and_check_INT(param_table_handle, "debug",
			   false, 0,	/* optional parameter defaulting to 0 */
			   false, 0, 0, NULL,
			   &debug);
if (status != 0)
   then return status;					/*** ERROR RETURN ***/

  {
CCTK_INT order;
status = get_and_check_INT(param_table_handle, "order",
			   true, 0,	/* this is a mandatory parameter */
			   true, 1, MAX_ORDER, "MAX_ORDER",   /* range check */
			   &order);
if (status != 0)
   then return status;					/*** ERROR RETURN ***/

  {
const int N_boundaries = 2*N_dims;
CCTK_INT  N_boundary_points_to_omit[MAX_N_BOUNDARIES];
status = get_INT_array(param_table_handle, "N_boundary_points_to_omit",
		       true, 0,		/* default value */
		       N_boundaries, N_boundary_points_to_omit,
		       NULL);
if (status != 0)
   then return status;					/*** ERROR RETURN ***/

/**************************************/

  {
CCTK_REAL boundary_off_centering_tolerance[MAX_N_BOUNDARIES];
CCTK_REAL boundary_extrapolation_tolerance[MAX_N_BOUNDARIES];
status = get_REAL_array(param_table_handle, "boundary_off_centering_tolerance",
			boundary_off_cntr_tol_default,
			N_boundaries, boundary_off_centering_tolerance);
if (status != 0)
   then return status;					/*** ERROR RETURN ***/
status = get_REAL_array(param_table_handle, "boundary_extrapolation_tolerance",
			boundary_extrap_tol_default,
			N_boundaries, boundary_extrapolation_tolerance);
if (status != 0)
   then return status;					/*** ERROR RETURN ***/
check_boundary_tolerances(N_boundaries,
			  boundary_off_centering_tolerance,
			  boundary_extrapolation_tolerance);

/**************************************/

  {
struct error_info error_info;
status = get_error_point_info(param_table_handle,
			      &error_info);
if (status != 0)
   then return status;					/*** ERROR RETURN ***/

/**************************************/

  {
#define MOLECULE_FAMILY_BUFSIZ	(MAX_MOLECULE_FAMILY_STRLEN+1)
char molecule_family_string[MOLECULE_FAMILY_BUFSIZ];
enum molecule_family molecule_family;
status = get_and_decode_molecule_family
		(param_table_handle,
		 MOLECULE_FAMILY_BUFSIZ, molecule_family_string,
		 &molecule_family);
if (status != 0)
   then return status;					/*** ERROR RETURN ***/

  {
CCTK_INT smoothing;
status = get_and_check_INT(param_table_handle, "smoothing",
			   false, 0,	/* optional, default value 0 */
			   true, 0, MAX_SMOOTHING, "MAX_SMOOTHING",
							/* range check */
			   &smoothing);
if (status != 0)
   then return status;					/*** ERROR RETURN ***/

/**************************************/

  {
CCTK_INT* const input_array_offsets
	= (CCTK_INT*) malloc(N_input_arrays1 * sizeof(CCTK_INT));
if (input_array_offsets == NULL)
   then return UTIL_ERROR_NO_MEMORY;			/*** ERROR RETURN ***/
status = get_INT_array(param_table_handle, "input_array_offsets",
		       true, 0,		/* default value */
		       N_input_arrays, input_array_offsets,
		       NULL);
if (status != 0)
   then {
	free(input_array_offsets);
	return status;					/*** ERROR RETURN ***/
	}

  {
CCTK_INT input_array_strides[MAX_N_DIMS];
status = get_INT_array(param_table_handle, "input_array_strides",
		       false, 0,	/* don't set default value */
		       N_dims, input_array_strides,
		       &value_not_set);
if (status != 0)
   then {
	free(input_array_offsets);
	return status;					/*** ERROR RETURN ***/
	}
if (value_not_set)
   then {
	/* default stride = Fortran storage order, */
	/* determined from input_array_dims[] */
	if (input_array_dims == NULL)
	   then {
		free(input_array_offsets);
		CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
			   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): can't default \"input_array_strides\"\n"
"                              with input_array_dims == NULL!");
		return UTIL_ERROR_BAD_INPUT;		/*** ERROR RETURN ***/
		}
	  {
	int stride = 1;
	int axis;
		for (axis = 0 ; axis < N_dims ; ++axis)
		{
		input_array_strides[axis] = stride;
		stride *= input_array_dims[axis];
		}
	  }
	}

  {
CCTK_INT input_array_min_subscripts[MAX_N_DIMS];
CCTK_INT input_array_max_subscripts[MAX_N_DIMS];
status = get_INT_array(param_table_handle, "input_array_min_subscripts",
		       true, 0,		/* default value */
		       N_dims, input_array_min_subscripts,
		       NULL);
if (status != 0)
   then {
	free(input_array_offsets);
	return status;					/*** ERROR RETURN ***/
	}

status = get_INT_array(param_table_handle, "input_array_max_subscripts",
		       false, 0,	/* don't set default value */
		       N_dims, input_array_max_subscripts,
		       &value_not_set);
if (status != 0)
   then {
	free(input_array_offsets);
	return status;					/*** ERROR RETURN ***/
	}
if (value_not_set)
   then {
	/* default max subscript = input_array_dims[]-1 along each axis */
	if (input_array_dims == NULL)
	   then {
		free(input_array_offsets);
		CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
			   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): can't default \"input_array_max_subscripts\"\n"
"                              with input_array_dims == NULL!");
		return UTIL_ERROR_BAD_INPUT;		/*** ERROR RETURN ***/
		}
	  {
	int axis;
		for (axis = 0 ; axis < N_dims ; ++axis)
		{
		input_array_max_subscripts[axis] = input_array_dims[axis] - 1;
		}
	  }
	}

/**************************************/

  {
CCTK_INT* const operand_indices
	= (CCTK_INT*) malloc(N_output_arrays1 * sizeof(CCTK_INT));
if (operand_indices == NULL)
   then {
	free(input_array_offsets);
	return UTIL_ERROR_NO_MEMORY;			/*** ERROR RETURN ***/
	}
status = get_INT_array(param_table_handle, "operand_indices",
		       false, 0,	/* don't set default value */
		       N_output_arrays, operand_indices,
		       &value_not_set);
if (status != 0)
   then {
	free(input_array_offsets);
	free(operand_indices);
	return UTIL_ERROR_NO_MEMORY;			/*** ERROR RETURN ***/
	}
if (value_not_set)
   then {
	/* default operand will use each input exactly once, */
	/* but this only makes since if N_input_arrays == N_output_arrays */
	if (N_input_arrays != N_output_arrays)
	   then {
		free(input_array_offsets);
		free(operand_indices);
		CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
			   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): can't default \"operand_indices\"\n"
"                              with N_input_arrays=%d != N_output_arrays=%d!"
			   ,
			   N_input_arrays, N_output_arrays);
		return UTIL_ERROR_BAD_INPUT;		/*** ERROR RETURN ***/
		}
	  {
	int out;
		for (out = 0 ; out < N_output_arrays ; ++out)
		{
		operand_indices[out] = out;
		}
	  }
	}

/**************************************/

  {
CCTK_INT* const operation_codes
	= (CCTK_INT*) malloc(N_output_arrays1 * sizeof(CCTK_INT));
if (operation_codes == NULL)
   then {
	free(input_array_offsets);
	free(operand_indices);
	return UTIL_ERROR_NO_MEMORY;			/*** ERROR RETURN ***/
	}
status = get_INT_array(param_table_handle, "operation_codes",
		       true, 0,		/* default value */
		       N_output_arrays, operation_codes,
		       NULL);
if (status != 0)
   then {
	free(input_array_offsets);
	free(operand_indices);
	free(operation_codes);
	return UTIL_ERROR_NO_MEMORY;			/*** ERROR RETURN ***/
	}

/******************************************************************************/

/*
 * set up for any molecule min/max m and/or position queries
 */

  {
struct molecule_min_max_m_info* p_molecule_min_max_m_info = NULL;
struct molecule_min_max_m_info  molecule_min_max_m_info;

/* molecule min/max m */
status1 = Util_TableQueryValueInfo(param_table_handle,
				   NULL, NULL,
				   "molecule_min_m");
status2 = Util_TableQueryValueInfo(param_table_handle,
				   NULL, NULL,
				   "molecule_max_m");
if ((status1 < 0) || (status2 < 0))
   then {
	free(input_array_offsets);
	free(operand_indices);
	free(operation_codes);
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): bad \"molecule_{min,max}_m\"\n"
"                              table entry/entries in query!\n"
"                              error status1=%d status2=%d"
		   ,
		   status1, status2);
	return (status1 < 0) ? status1 : status2;	/*** ERROR RETURN ***/
	}
if (status1 && status2)
   then p_molecule_min_max_m_info = &molecule_min_max_m_info;

  {
CCTK_INT** p_molecule_positions = NULL;
CCTK_INT*  molecule_positions_array[MAX_N_DIMS];

/* are we doing a molecule-positions query? */
status = Util_TableQueryValueInfo(param_table_handle,
				  NULL, NULL,
				  "molecule_positions");
if (status < 0)
   then {
	free(input_array_offsets);
	free(operand_indices);
	free(operation_codes);
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): bad \"molecule_positions\"\n"
"                              table entry in query!\n"
"                              error status=%d"
		   ,
		   status);
	return status;					/*** ERROR RETURN ***/
	}
if (status)
   then {
	/* yes, we're doing a molecule-positions query */
	status = get_molecule_positions(param_table_handle,
				       N_dims, molecule_positions_array);
	if (status != 0)
	   then return status;				/*** ERROR RETURN ***/
	p_molecule_positions = molecule_positions_array;
	}

/******************************************************************************/

/*
 * set up for any Jacobian queries
 */

  {
struct Jacobian_info* p_Jacobian_info = NULL;
struct Jacobian_info  Jacobian_info;
Jacobian_info.Jacobian_pointer = NULL;
Jacobian_info.Jacobian_offset = NULL;

/* are we doing a Jacobian query? */
status = Util_TableQueryValueInfo(param_table_handle,
				  NULL, NULL,
				  "Jacobian_pointer");
if (status < 0)
   then {
	free(input_array_offsets);
	free(operand_indices);
	free(operation_codes);
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): bad \"Jacobian_pointer\" table entry!\n"
"                              (preliminary check)\n"
"                              error status=%d"
		   ,
		   status);
	 return status;					/*** ERROR RETURN ***/
	}

if (status)
   then {
	/* yes, we're doing a Jacobian query */
	status = get_Jacobian_info(param_table_handle,
				   N_dims, N_output_arrays,
				   &Jacobian_info);
	if (status != 0)
	   then {
		free(input_array_offsets);
		free(operand_indices);
		free(operation_codes);
		return status;				/*** ERROR RETURN ***/
		}
	p_Jacobian_info = & Jacobian_info;
	}

/******************************************************************************/

/*
 * decode (interp_operator, N_dims, molecule_family, order, smoothing)
 * and call the appropriate subfunction to do the actual interpolation
 */


/* firewall: make sure all the lookup indices are in range */
/* ... commented-out comparisons are always true since enums are unsigned */
if (    /*(interp_operator >= 0) &&*/ (interp_operator <  N_INTERP_OPERATORS)
     &&   (N_dims          >= 0) &&   (N_dims          <= MAX_N_DIMS)
     && /*(molecule_family >= 0) &&*/ (molecule_family <  N_MOLECULE_FAMILIES)
     &&   (order           >= 0) &&   (order           <= MAX_ORDER)
     &&   (smoothing       >= 0) &&   (smoothing       <= MAX_SMOOTHING)    )
   then {
	/* ok ==> no-op here */
	}
   else {
	CCTK_CVWarn(BUG_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform():\n"
"        internal error (interpolator bug!):\n"
"        one or more function-pointer-table indices is out of range!\n"
"        interp_operator=(int)%d, should be in [0,%d)\n"
"        N_dims         =     %d, should be in [0,%d]\n"
"        molecule_family=(int)%d, should be in [0,%d)\n"
"        order          =     %d, should be in [0,%d]\n"
"        smoothing      =     %d, should be in [0,%d]\n"
			   ,
			   (int) interp_operator, (int) N_INTERP_OPERATORS,
			   (int) N_dims         , (int) MAX_N_DIMS,
			   (int) molecule_family, (int) N_MOLECULE_FAMILIES,
			   (int) order          , (int) MAX_ORDER,
			   (int) smoothing      , (int) MAX_SMOOTHING);
								/*NOTREACHED*/
	}

/* look up the subfunction to do the interpolation */
  {
const p_interp_fn_t p_interp_fn = p_interp_fn_table[interp_operator]
						   [N_dims]
						   [molecule_family]
						   [order]
						   [smoothing];
if (p_interp_fn == NULL_INTERP_FN_PTR)
   then {
	free(input_array_offsets);
	free(operand_indices);
	free(operation_codes);
	if (p_Jacobian_info != NULL)
	   then {
		free(Jacobian_info.Jacobian_offset);
		free(Jacobian_info.Jacobian_pointer);
		}
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform():\n"
"        interpolation not implemented for the combination\n"
"        interp_operator=\"%s\", N_dims=%d\n"
"        molecule_family=\"%s\", order=%d, smoothing=%d",
		   interp_operator_string, N_dims,
		   molecule_family_string, (int)order, (int)smoothing);
	return UTIL_ERROR_BAD_INPUT;			/*** ERROR RETURN ***/
	}

/* call the subfunction to actually do the interpolation */
if (debug > 0)
   then {
	printf("AEILocalInterp::InterpLocalUniform.c: calling interpolator fn (N_interp_points=%d)\n", N_interp_points);
	fflush(stdout);
	}
  {
struct molecule_structure_flags molecule_structure_flags;
const int return_code
	= (*p_interp_fn)(coord_origin, coord_delta,
			 N_interp_points,
			    interp_coords_type_code, interp_coords,
			    N_boundary_points_to_omit,
			    boundary_off_centering_tolerance,
			    boundary_extrapolation_tolerance,
			 N_input_arrays,
			    input_array_offsets, input_array_strides,
			    input_array_min_subscripts,
			    input_array_max_subscripts,
			    input_array_type_codes, input_arrays,
			 N_output_arrays,
			    output_array_type_codes, output_arrays,
			    operand_indices, operation_codes,
			 debug, (false ? NULL : NULL),
			 &error_info,
			 &molecule_structure_flags,
			 p_molecule_min_max_m_info,
			 p_molecule_positions,
			 p_Jacobian_info);
if (debug > 0)
   then {
	printf("AEILocalInterp::InterpLocalUniform.c: back from interpolator fn with return_code=%d\n", return_code);
	fflush(stdout);
	}

/******************************************************************************/

/*
 * set any further error status in the parameter table
 */
status = set_error_info(param_table_handle,
			&error_info);
if (status != 0)
   then {
	free(input_array_offsets);
	free(operand_indices);
	free(operation_codes);
	if (p_Jacobian_info != NULL)
	   then {
		free(Jacobian_info.Jacobian_offset);
		free(Jacobian_info.Jacobian_pointer);
		}
	return status;					/*** ERROR RETURN ***/
	}

/******************************************************************************/

/*
 * store query results
 *
 * ... only molecule structure and min/max m have to be stored explicitly;
 *     the other queries are stored directly into their final destinatios
 *     by the interpolation subfunction
 */

status = set_molecule_structure(param_table_handle,
				&molecule_structure_flags);
if (status != 0)
   then {
	free(input_array_offsets);
	free(operand_indices);
	free(operation_codes);
	if (p_Jacobian_info != NULL)
	   then {
		free(Jacobian_info.Jacobian_offset);
		free(Jacobian_info.Jacobian_pointer);
		}
	return status;					/*** ERROR RETURN ***/
	}

if (p_molecule_min_max_m_info != NULL)
   then {
	status = set_molecule_min_max_m(param_table_handle,
					N_dims, p_molecule_min_max_m_info);
	if (status != 0)
	   then {
		free(input_array_offsets);
		free(operand_indices);
		free(operation_codes);
		if (p_Jacobian_info != NULL)
		   then {
			free(Jacobian_info.Jacobian_offset);
			free(Jacobian_info.Jacobian_pointer);
			}
		return status;				/*** ERROR RETURN ***/
		}
	}

/******************************************************************************/

free(input_array_offsets);
free(operand_indices);
free(operation_codes);
if (p_Jacobian_info != NULL)
   then {
	free(Jacobian_info.Jacobian_offset);
	free(Jacobian_info.Jacobian_pointer);
	}

if (debug > 0)
   then {
	printf("AEILocalInterp::InterpLocalUniform.c: returning %d\n", return_code);
	fflush(stdout);
	}
return return_code;					/*** NORMAL RETURN ***/
  }
  }
  }
  }
  }
  }
  }
  }
  }
  }
  }
  }
  }
  }
  }
  }
  }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*
 * This function checks
 *	CCTK_REAL boundary_off_centering_tolerance[N_boundaries]
 *	CCTK_REAL boundary_extrapolation_tolerance[N_boundaries]
 * for validity, warning about any dubious settings.
 */
static
  void check_boundary_tolerances
	(int N_boundaries,
	 const CCTK_REAL boundary_off_centering_tolerance[MAX_N_BOUNDARIES],
	 const CCTK_REAL boundary_extrapolation_tolerance[MAX_N_BOUNDARIES])
{
int ibndry;
	for (ibndry = 0 ; ibndry < N_boundaries ; ++ibndry)
	{
	if (    (boundary_off_centering_tolerance[ibndry] == 0.0)
	     && (boundary_extrapolation_tolerance[ibndry] > 0.0)    )
	   then CCTK_CVWarn(WARNING_MSG_SEVERITY_LEVEL,
			   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform():\n"
"        warning: setting boundary_off_centering_tolerance[%d] == 0.0\n"
"                 and     boundary_extrapolation_tolerance[%d] == %g > 0.0\n"
"                 is almost certainly a mistake\n"
"                 (the boundary_extrapolation_tolerance[%d] == %g\n"
"                  setting will have no effect because of the\n"
"                  boundary_off_centering_tolerance[%d] == 0.0 setting)"
			   ,
			   ibndry,
			   ibndry,
			      (double) boundary_extrapolation_tolerance[ibndry],
			   ibndry,
			      (double) boundary_extrapolation_tolerance[ibndry],
			   ibndry);
	}
}

/******************************************************************************/

/*
 * This function tries to get
 *	CCTK_POINTER per_point_status
 * from the parameter table, and sets  found_per_point_status  to
 * indicate whether or not that key was found.
 *
 * If  per_point_status  is found, this function casts it to a CCTK_INT*
 * and stores it in p_error_info->per_point_status.  If not (if there's
 * no such key in the table), this function sets
 * p_error_info->per_point_status to a NULL pointer.
 *
 * Results:
 * This function returns 0 for ok, or the (nonzero) return code for error.
 */
static
  int get_error_point_info(int param_table_handle,
			   struct error_info *p_error_info)
{
CCTK_POINTER per_point_status;
int status;

status = Util_TableGetPointer(param_table_handle,
			      &per_point_status,
			      "per_point_status");
if      (status == 1)
   then {
	p_error_info->found_per_point_status = true;
	p_error_info->per_point_status = (CCTK_INT*) per_point_status;
	}
else if (status == UTIL_ERROR_TABLE_NO_SUCH_KEY)
   then {
	p_error_info->found_per_point_status = false;
	p_error_info->per_point_status = NULL;
	}
else	{
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): bad \"per_point_status\" table entry!\n"
"                              error status=%d",
		   status);
	return status;					/*** ERROR RETURN ***/
	}

status = Util_TableQueryValueInfo(param_table_handle,
				  NULL, NULL,
				  "suppress_warnings");
if	(status == 1)
   then {
	/* key "suppress_warnings" is in table ==> no warning msgs */
	p_error_info->print_warning_msg = false;
	}
else if (status == 0)
   then {
	/* key "suppress_warnings" is NOT in table */
	/* ==> leave warning msgs on by default*/
	p_error_info->print_warning_msg = true;
	}
else	{
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): bad \"suppress_warnings\" table entry!\n"
"                              error status=%d",
		   status);
	return status;					/*** ERROR RETURN ***/
	}

return 0;						/*** NORMAL RETURN ***/
}

/******************************************************************************/

/*
 * This function gets
 *	const char molecule_family[]
 * from the parameter table, or sets the default value "cube" if this
 * key isn't present.  This function also decodes  molecule_family  into
 * an  enum molecule_family .
 *
 * Results:
 * This function returns 0 for ok, or the (nonzero) return code for error.
 */
static
  int get_and_decode_molecule_family
	(int param_table_handle,
	 int buffer_size, char molecule_family_string_buffer[],
	 enum molecule_family *p_molecule_family)
{
enum molecule_family molecule_family;
int status;

status = Util_TableGetString(param_table_handle,
			     buffer_size, molecule_family_string_buffer,
			     "molecule_family");

if	(status == UTIL_ERROR_TABLE_NO_SUCH_KEY)
   then {
	/* set the default value */
	Util_Strlcpy(molecule_family_string_buffer, "cube", buffer_size);
	molecule_family = molecule_family_cube;

	/* set this key in the parameter table */
	/* to give the (default) molecule family we're going to use */
	status = Util_TableSetString(param_table_handle,
				     molecule_family_string_buffer,
				     "molecule_family");
	if (status < 0)
	   then {
		CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
			   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): error setting default molecule family\n"
"                              in parameter table!\n"
"                              error status=%d",
			   status);
		return status;				/*** ERROR RETURN ***/
		}
	}
else if (status > 0)
   then {
	/* decode molecule family string */
	if (strcmp(molecule_family_string_buffer, "cube") == 0)
	   then molecule_family = molecule_family_cube;
	else	{
		CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
			   __LINE__, __FILE__, CCTK_THORNSTRING,
"CCTK_InterpLocalUniform(): unknown molecule_family string \"%s\"!",
			   molecule_family_string_buffer);
		return UTIL_ERROR_BAD_INPUT;		/*** ERROR RETURN ***/
		}
	}
else	{
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): bad table entry for\n"
"                              \"molecule_family\" parameter!\n"
"                              error status=%d",
		   status);
	return status;					/*** ERROR RETURN ***/
	}

*p_molecule_family = molecule_family;
return 0;						/*** NORMAL RETURN ***/
}

/******************************************************************************/

/*
 * This function gets
 *	CCTK_INT* molecule_positions[N_dims]
 * from the parameter table.
 *
 * Results:
 * This function returns 0 for ok, or the (nonzero) return code for error.
 */
static
  int get_molecule_positions
	(int param_table_handle,
	 int N_dims, CCTK_INT* molecule_positions_array[MAX_N_DIMS])
{
CCTK_POINTER molecule_positions_temp[MAX_N_DIMS];
int status;

status = Util_TableGetPointerArray(param_table_handle,
				   N_dims, molecule_positions_temp,
				   "molecule_positions");

if (status == N_dims)
   then {
	/* type-cast CCTK_POINTER pointers to CCTK_INT* pointers */
	/* (which point to where the query results will be stored) */
	int axis;
		for (axis = 0 ; axis < N_dims ; ++axis)
		{
		molecule_positions_array[axis]
			= (CCTK_INT *) molecule_positions_temp[axis];
		}
	}
else	{
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): bad \"molecule_positions\" table entry!\n"
"                              error status=%d",
		   status);
	return status;					/*** ERROR RETURN ***/
	}

return 0;						/*** NORMAL RETURN ***/
}

/******************************************************************************/

/*
 * This function gets the Jacobian-query parameters
 *	CCTK_REAL* Jacobian_pointer[N_output_arrays]
 *	CCTK_INT   Jacobian_offset [N_output_arrays]   # optional, default=all 0
 *	CCTK_INT   Jacobian_interp_point_stride
 *	CCTK_INT   Jacobian_m_point_strides[N_dims]
 *	CCTK_INT   Jacobian_part_stride                # optional, default=1
 * and stores them in a  struct Jacobian_info .
 *
 * Results:
 * This function returns 0 for ok, or the (nonzero) return code for error.
 */
static
  int get_Jacobian_info(int param_table_handle,
			int N_dims, int N_output_arrays,
			struct Jacobian_info* p_Jacobian_info)
{
/* padded array size, cf. InterpLocalUniform() header comments */
const int N_output_arrays1 = N_output_arrays + 1;

int status;

p_Jacobian_info->Jacobian_pointer
	= (CCTK_REAL**) malloc(N_output_arrays1 * sizeof(CCTK_REAL *));
if (p_Jacobian_info->Jacobian_pointer == NULL)
   then return UTIL_ERROR_NO_MEMORY;			/*** ERROR RETURN ***/

  {
CCTK_POINTER* Jacobian_pointer_temp
	= (CCTK_POINTER*) malloc(N_output_arrays1 * sizeof(CCTK_POINTER));
if (Jacobian_pointer_temp == NULL)
   then {
	free(p_Jacobian_info->Jacobian_pointer);
	return UTIL_ERROR_NO_MEMORY;			/*** ERROR RETURN ***/
	}
status = Util_TableGetPointerArray(param_table_handle,
				   N_output_arrays, Jacobian_pointer_temp,
				   "Jacobian_pointer");
if	(status == N_output_arrays)
   then {
	/* type-cast CCTK_POINTER pointers to CCTK_REAL* pointers */
	/* (which point to where the query results will be stored) */
	int out;
		for (out = 0 ; out < N_output_arrays ; ++out)
		{
		p_Jacobian_info->Jacobian_pointer[out]
			= (CCTK_REAL *) Jacobian_pointer_temp[out];
		}
	free(Jacobian_pointer_temp);
	}
else	{
	free(p_Jacobian_info->Jacobian_pointer);
	free(Jacobian_pointer_temp);
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): bad \"Jacobian_pointer\" table entry!\n"
"                              error status=%d"
		   ,
		   status);
	return status;					/*** ERROR RETURN ***/
	}
  }

p_Jacobian_info->Jacobian_offset
	= (CCTK_INT*) malloc(N_output_arrays1 * sizeof(CCTK_INT));
if (p_Jacobian_info->Jacobian_offset == NULL)
   then {
	free(p_Jacobian_info->Jacobian_pointer);
	return UTIL_ERROR_NO_MEMORY;			/*** ERROR RETURN ***/
	}
status = get_INT_array(param_table_handle, "Jacobian_offset",
		       true, 0,		/* default value */
		       N_output_arrays, p_Jacobian_info->Jacobian_offset,
		       NULL);
if (status != 0)
   then {
	free(p_Jacobian_info->Jacobian_pointer);
	free(p_Jacobian_info->Jacobian_offset);
	return status;					/*** ERROR RETURN ***/
	}

status = get_and_check_INT(param_table_handle, "Jacobian_interp_point_stride",
			   true, 0,	/* this is a mandatory parameter */
					/* if we're executing this function */
			   false, 0, 0, NULL,	/* no range check */
			   & p_Jacobian_info->Jacobian_interp_point_stride);
if (status != 0)
   then {
	free(p_Jacobian_info->Jacobian_offset);
	free(p_Jacobian_info->Jacobian_pointer);
	return status;					/*** ERROR RETURN ***/
	}

status = get_INT_array(param_table_handle, "Jacobian_m_strides",
		       false, 0,	/* don't set default value */
		       N_dims, p_Jacobian_info->Jacobian_m_strides,
		       NULL);		/* this is a mandatory parameter */
					/* if we're executing this function */
if (status != 0)
   then {
	free(p_Jacobian_info->Jacobian_pointer);
	free(p_Jacobian_info->Jacobian_offset);
	return status;					 /*** ERROR RETURN ***/
	}

status = get_and_check_INT(param_table_handle, "Jacobian_part_stride",
			   false, 1,	/* default value */
			   false, 0, 0, NULL,	/* no range check */
			   & p_Jacobian_info->Jacobian_part_stride);
if (status != 0)
   then {
	free(p_Jacobian_info->Jacobian_pointer);
	free(p_Jacobian_info->Jacobian_offset);
	return status;					/*** ERROR RETURN ***/
	}

return 0;						/*** NORMAL RETURN ***/
}

/******************************************************************************/

/*
 * If p_error_info->found_per_point_status is true, this function
 * sets the paramater-table entry
 *	CCTK_INT error_point_status
 * to the negative of p_error_info->error_count.
 *
 * Results:
 * This function returns 0 for ok, or the (nonzero) return code for error.
 */
static
  int set_error_info(int param_table_handle,
		     struct error_info* p_error_info)
{
if (p_error_info->found_per_point_status)
   then {
	const int status = Util_TableSetInt(param_table_handle,
					    -p_error_info->error_count,
					    "error_point_status");
	if (status < 0)
	   then {
		CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
			   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform():\n"
"        error setting \"error_point_status\" in parameter table!"
"        error status=%d"
			   ,
			   status);				
		return status;				/*** ERROR RETURN ***/
		}
	}

return 0;						/*** NORMAL RETURN ***/
}

/******************************************************************************/

/*
 * This function sets the molecule structure flags
 *	CCTK_INT MSS_is_fn_of_interp_coords
 *	CCTK_INT MSS_is_fn_of_which_operation
 *	CCTK_INT MSS_is_fn_of_input_array_values
 *	CCTK_INT Jacobian_is_fn_of_input_array_values
 * in the parameter table.
 *
 * Results:
 * This function returns 0 for ok, or the (nonzero) return code for error.
 */
static
  int set_molecule_structure
	(int param_table_handle,
	 const struct molecule_structure_flags* p_molecule_structure_flags)
{
const int status1
   = Util_TableSetInt(param_table_handle,
		      p_molecule_structure_flags->MSS_is_fn_of_interp_coords,
		      "MSS_is_fn_of_interp_coords");
const int status2
   = Util_TableSetInt(param_table_handle,
		      p_molecule_structure_flags->MSS_is_fn_of_which_operation,
		      "MSS_is_fn_of_which_operation");
const int status3
   = Util_TableSetInt(param_table_handle,
		      p_molecule_structure_flags->MSS_is_fn_of_input_array_values,
		      "MSS_is_fn_of_input_array_values");
const int status4
   = Util_TableSetInt(param_table_handle,
		      p_molecule_structure_flags->Jacobian_is_fn_of_input_array_values,
		      "Jacobian_is_fn_of_input_array_values");
if ((status1 < 0) || (status2 < 0) || (status3 < 0) || (status4 < 0))
   then {
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform():\n"
"        error setting molecule structure flags table entry/entries!"
"        error status1=%d status2=%d status3=%d status4=%d"
			   ,
			   status1, status2, status3, status4);
	return   (status1 < 0) ? status1		/*** ERROR RETURN ***/
	       : (status2 < 0) ? status2		/*** ERROR RETURN ***/
	       : (status3 < 0) ? status3		/*** ERROR RETURN ***/
	       :		 status4;		/*** ERROR RETURN ***/
	}

return 0;						/*** NORMAL RETURN ***/
}

/******************************************************************************/

/*
 * This function sets the molecule size parameters
 *	CCTK_INT molecule_min_m[N_dims]
 *	CCTK_INT molecule_max_m[N_dims]
 * in the parameter table.
 *
 * Results:
 * This function returns 0 for ok, or the (nonzero) return code for error.
 */
static
  int set_molecule_min_max_m
	(int param_table_handle,
	 int N_dims,
	 const struct molecule_min_max_m_info* p_molecule_min_max_m_info)
{
const int status1
   = Util_TableSetIntArray(param_table_handle,
			   N_dims, p_molecule_min_max_m_info->molecule_min_m,
			   "molecule_min_m");
const int status2
   = Util_TableSetIntArray(param_table_handle,
			   N_dims, p_molecule_min_max_m_info->molecule_max_m,
			   "molecule_max_m");
if ((status1 < 0) || (status2 < 0))
   then {
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): error setting\n"
"                              \"molecule_{min,max}_m\" table entry/entries!"
"                              error status1=%d status2=%d",
		   status1, status2);
	return (status1 < 0) ? status1 : status2;	/*** ERROR RETURN ***/
	}

return 0;						/*** NORMAL RETURN ***/
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*
 * This function gets and range-checks a CCTK_INT from the parameter table.
 *
 * Arguments:
 * param_table_handle = handle to the Cactus key-value table
 * name = the character-string key in the table
 * mandatory_flag = true  ==> the value must be present in the parameter
 *			      table (default_value is ignored)
 *		    false ==> the value is optional
 * default_value = the default value (for each array element) if the
 *		   key isn't in the parameter table and mandatory_flag == false
 * check_range_flag = true  ==> check that the value satisfies
 *				min_value <= value <= max_value
 *		      false ==> don't do a range check, and ignore
 *				min_value, max_value, and max_value_string
 * {min,max}_value = the inclusive range of valid values
 * max_value_string = the character-string name of max_value (this is used
 *		       only in formatting the out-of-range error message)
 * p_value --> where we should store the value
 *
 * Results:
 * This function returns 0 for ok, or the (nonzero) return code for error.
 */
static
  int get_and_check_INT(int param_table_handle, const char name[],
			bool mandatory_flag, int default_value,
			bool check_range_flag, int min_value, int max_value,
			const char max_value_string[],
			CCTK_INT* p_value)
{
CCTK_INT value;

const int status = Util_TableGetInt(param_table_handle, &value, name);

if	(status == 1)
   then { }	/* value set from table ==> no-op here */
else if (status == UTIL_ERROR_TABLE_NO_SUCH_KEY)
   then {
	if (mandatory_flag)
	   then {
		CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
			   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): missing table entry for\n"
"                              \"%s\" parameter!\n"
"                              (this is a mandatory parameter)!\n"
			   ,
			   name);
		return status;				/*** ERROR RETURN ***/
		}
	   else value = default_value;
	}
else	{
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): bad table entry for\n"
"                              \"%s\" parameter!\n"
"                              error status=%d"
		   ,
		   name,
		   status);
	return status;					/*** ERROR RETURN ***/
	}

if (check_range_flag)
   then {
	if ((value < min_value) || (value > max_value))
	   then {
		CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
			   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): %s=%d is invalid!\n"
"                              valid range is %d <= %s <= %s=%d",
			   name, (int)value,
			   min_value, name, max_value_string, max_value);
		return UTIL_ERROR_BAD_INPUT;		/*** ERROR RETURN ***/
		}
	}

*p_value = value;
return 0;						/*** NORMAL RETURN ***/
}

/******************************************************************************/

/*
 * This function gets an array of CCTK_INT values from the parameter table.
 *
 * Arguments:
 * param_table_handle = handle to the Cactus key-value table
 * name = the character-string key in the table
 * default_flag = false ==> ignore  default_value  and set
 *			     *p_value_not_set  to true if the
 *			    key is not present in the parameter table,
 *			    false if it is present
 *		  true  ==> if the key is not present in the parameter table,
 *			    set each array element to  default_value
 * default_value = the default value (for each array element) if the
 *		   key isn't in the parameter table
 * N = the size of the array
 * buffer --> a buffer in which to store the array
 * p_value_not_set = NULL or a pointer to a flag to be set if and only if
 *		     this function did *not* set values in the buffer.
 *		     More precisely, if  default_flag  is set, this argument
 *		     is ignored.  Otherwise, the logic is this:
 *			if (this function set values in the buffer)
 *			   then {
 *				if (p_value_not_set is non-NULL)
 *				   then *p_value_not_set = false
 *				}
 *			   else {
 *				if (p_value_not_set is non-NULL)
 *				   then *p_value_not_set = true
 *				   else give an error message that
 *					a mandatory value is missing
 *					from the parameter table, then
 *					return UTIL_ERROR_TABLE_NO_SUCH_KEY
 *				}
 *
 * Results:
 * This function returns 0 for ok, or the desired (nonzero) error return
 * code if something is wrong.  Note that error return codes may be either
 * positive or negative!
 */
static
  int get_INT_array(int param_table_handle, const char name[],
		    bool default_flag, int default_value,
		    int N, CCTK_INT buffer[],
		    bool* p_value_not_set)
{
const int status
	= Util_TableGetIntArray(param_table_handle,
				N, buffer,
				name);

if	 (status == UTIL_ERROR_TABLE_NO_SUCH_KEY)
   then {
	/* we didn't set values in the buffer */
	/* ==> either set the default value ourself, */
	/*     or print an error message, */
	/*     or notify our caller that we didn't do so */
	if (default_flag)
	   then {
		int i;
			for (i = 0 ; i < N ; ++i)
			{
			buffer[i] = default_value;
			}
		}
	   else {
		if (p_value_not_set == NULL)
		   then {
			CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
				   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): missing table entry for\n"
"                              \"%s\" parameter!\n"
"                              (this is a mandatory parameter)!\n"
				   ,
				   name);
			return UTIL_ERROR_TABLE_NO_SUCH_KEY;  /*** ERROR RETURN ***/
			}
		   else *p_value_not_set = true;
		}
	}
else if (status == N)
   then {
	/* we did set values in the buffer */
	if (! default_flag)
	   then if (p_value_not_set != NULL)
		   then *p_value_not_set = false;
	}
else	{
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): bad table entry for\n"
"                              \"%s\" parameter!\n"
"                              error status=%d",
		   name, status);
	return status;					/*** ERROR RETURN ***/
	}

return 0;						/*** NORMAL RETURN ***/
}

/******************************************************************************/

/*
 * This function gets an array of CCTK_REAL values from the parameter table.
 *
 * Arguments:
 * param_table_handle = handle to the Cactus key-value table
 * name = the character-string key in the table
 * default_value = the default value (for each array element) if the
 *		   key isn't in the parameter table
 * N = the size of the array
 * buffer --> a buffer in which to store the array
 *
 * Results:
 * This function returns 0 for ok, or the desired (nonzero) error return
 * code if something is wrong.  Note that error return codes may be either
 * positive or negative!
 */
static
  int get_REAL_array(int param_table_handle, const char name[],
		     CCTK_REAL default_value,
		     int N, CCTK_REAL buffer[])
{
const int status
	= Util_TableGetRealArray(param_table_handle,
				 N, buffer,
				 name);
if	 (status == UTIL_ERROR_TABLE_NO_SUCH_KEY)
   then {
	/* set default value */
	int i;
		for (i = 0 ; i < N ; ++i)
		{
		buffer[i] = default_value;
		}
	}
else if (status == N)
   then { }	/* ok ==> no-op here */
else	{
	CCTK_CVWarn(ERROR_MSG_SEVERITY_LEVEL,
		   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): bad table entry for\n"
"                              \"%s\" parameter!\n"
"                              error status=%d",
		   name, status);
	return status;					/*** ERROR RETURN ***/
	}

return 0;						/*** NORMAL RETURN ***/
}
