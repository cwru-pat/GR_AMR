/*@@
  @file      template.c
  @date      23.Oct.2001
  @author    Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc
	This file is a template for generalized interpolation for
	1-, 2-, or 3-d molecules.  It's used by defining various
	C preprocessor macros, then #including this file.  Each
	such inclusion defines a single interpolation function
	(plus other supporting functions local to this file),
	which does generalized interpolation for a single molecule,
	or more accurately for a single molecule for each possible
	operation_codes[] value.
  @enddesc

  @hdate    29.Aug.2002
  @hauthor  Jonathan Thornburg <jthorn@aei.mpg.de>
  @hdesc    Restructure code and break it up into smaller functions
	    to reduce excessive compilation cpu/memory with big
	    3-D molecules.
  @endhdesc

  @hdate    28.Jan.2003
  @hauthor  Jonathan Thornburg <jthorn@aei.mpg.de>
  @hdesc    Support parameter-table entries
            @var boundary_off_centering_tolerance and
            @var boundary_extrapolation_tolerance.
  @endhdesc 

  @desc
	The following header files must be #included before
	#including this file:
		<math.h>
		<limits.h>
		<stdlib.h>
		<string.h>
		<stdio.h>
		"util_ErrorCodes.h"
		"cctk.h"
		"InterpLocalUniform.h"

	The following preprocessor macros must be defined before
	#including this file (note the macros which are file names
	must all include the proper quotes for a #include in their
	expansion, e.g.
		#define DATA_VAR_DCL_FILE_NAME	\
			"coeffs/2d.cube.size3/data-var.dcl.c"
  @enddesc


  @var	    FUNCTION_NAME
  @vdesc    The name of the interpolation function, e.g.
		#define FUNCTION_NAME	AEILocalInterp_U_LagTP_2cube_20
  @endvar

  @var	    N_DIMS
  @vdesc    The number of dimensions in which to interpolate, e.g.
		#define N_DIMS	2
	    The present implementation restricts this to 1, 2,
	    or 3, but this could easily be changed if needed.
	    Note that MAX_N_DIMS (defined in "InterpLocalUniform.c")
	    is a compile-time upper bound for N_DIMS, useful for
	    sizing arrays etc.
  @endvar

  @var	    XYZ
  @desc	    A comma-separated list of the (x,y,z) coordinates in the
	    interpolation, e.g.
		#define XYZ	x,y
  @endvar

  @var	    FP_XYZ
  @desc	    A comma-separated list of function-prototype declarations
	    of the (x,y,z) coordinates in the interpolation, e.g.
		#define FP_XYZ	fp x, fp y
  @endvar

  @var	    STRIDE_IJK
  @desc	    A comma-separated list of the (i,j,k) array strides for the
	    interpolation, e.g.
		#define STRIDE_IJK	stride_i, stride_j
  @endvar

  @var	    JACOBIAN_MIJK_STRIDE
  @desc	    A comma-separated list of the (i,j,k) strides for the
	    Jacobian, e.g.
		#define JACOBIAN_MIJK_STRIDE	\
			Jacobian_mi_stride, Jacobian_mj_stride
  @endvar

  @var	    MOLECULE_MIN_M
  @vdesc    The minimum m coordinate in the molecule, e.g.
		#define MOLECULE_MIN_M	-1
  	    The present implementation takes this to be the same in each
	    dimension, but this could easily be changed if needed.  
  @endvar

  @var	    MOLECULE_MAX_M
  @vdesc    The maximum m coordinate in the molecule, e.g.
		#define MOLECULE_MAX_M	1
  	    The present implementation takes this to be the same in each
	    dimension, but this could easily be changed if needed.  
  @endvar

  @var	    MOLECULE_SIZE
  @vdesc    The diameter of (number of points in) the molecules to be used,
	    e.g.
		#define MOLECULE_SIZE	3
  	    The present implementation takes this to be the same in each
	    dimension, but this could easily be changed if needed.  
  @endvar

  @var	    HAVE_OP_{I,DX,DY,DXX,DXY,DYY,...}
  @vdesc    Each of these symbols should be defined or not defined
	    according as if the corresponding derivative operator is
	    to be supported or not supported by this function, e.g.
		#define HAVE_OP_I
		#define HAVE_OP_DX
		#define HAVE_OP_DY
		#define HAVE_OP_DXX
		#define HAVE_OP_DXY
		#define HAVE_OP_DYY
	    if we support all 1st and 2nd derivatives in 2-D.
  @endvar

  @var	    DATA_STRUCT
  @vdesc    The name of a C struct containing all the data variables,
	    e.g.
		#define DATA_STRUCT   data_struct_2d_cube_size3
  @endvar

  @var	    COEFF_STRUCT
  @vdesc    The name of a C struct containing all the molecule coefficients,
	    e.g.
		#define COEFFS_STRUCT coeffs_struct_2d_cube_size3
  @endvar

  @var	    LOAD_DATA_{REAL,REAL{4,8,16},COMPLEX,COMPLEX{8,16,32}}
  @vdesc    The name of a C function to load data from the input arrays
	    into the DATA_STRUCT structure.  Typically these will be
	    chosem from among the functions defined in
	    "common/load-template.[ch]", e.g.
		#define LOAD_DATA_REAL		AEILocalInterp_load_2dcube3_r
		#define LOAD_DATA_REAL4		AEILocalInterp_load_2dcube3_r4
		#define LOAD_DATA_REAL8		AEILocalInterp_load_2dcube3_r8
		#define LOAD_DATA_REAL16	AEILocalInterp_load_2dcube3_r16
		#define LOAD_DATA_COMPLEX	AEILocalInterp_load_2dcube3_c
		#define LOAD_DATA_COMPLEX8	AEILocalInterp_load_2dcube3_c8
		#define LOAD_DATA_COMPLEX16	AEILocalInterp_load_2dcube3_c16
		#define LOAD_DATA_COMPLEX32	AEILocalInterp_load_2dcube3_c32
  @endvar

  @var	    EVALUATE_MOLECULE
  @vdesc    The name of a C function to evaluate the dot product
	    of the DATA_STRUCT and COEFF_STRUCT structures.  Typically
	    these will be chosem from among the functions defined in
	    "common/evaluate.[ch]", e.g.
		#define EVALUATE_MOLECULE	LocalInterp_eval_2d_cube6
  @endvar

  @var	    STORE_COEFFS
  @vdesc    The name of a C function to store a COEDFF_STRUCT structure
	    of the DATA_STRUCT and COEFF_STRUCT structures.  Typically
	    this will be chosem from among the functions defined in
	    "common/store.[ch]", e.g.
		#define STORE_COEFFS		LocalInterp_store_2d_cube6
  @endvar

  @var	    COEFFS_{I,DX,DY,DXX,DXY,DYY,...}_COMPUTE_FILE_NAME
  @vdesc    Each of these macros should be the name of a file
	    (presumably machine-generated) containing a sequence
	    of C assignment statements to compute all the coefficients
	    in a corresponding  coeffs_*  structure as polynomials
	    in the variables (x,y,z), e.g. (for 2D size-3 molecules,
	    dx operator)
		fp t36;
		fp t42;
		fp t41;
		fp t40;
		fp t39;
		fp t37;
		      t36 = RATIONAL(1.0,3.0)*x;
		      t42 = RATIONAL(1.0,4.0)*y+t36;
		      t41 = t36+RATIONAL(-1.0,4.0)*y;
		      t40 = RATIONAL(1.0,6.0);
		      t39 = RATIONAL(-1.0,6.0);
		      t37 = RATIONAL(-2.0,3.0)*x;
		      coeffs_dx->coeff_m1_m1 = t39+t42;
		      coeffs_dx->coeff_0_m1 = t37;
		      coeffs_dx->coeff_p1_m1 = t40+t41;
		      coeffs_dx->coeff_m1_0 = t36+t39;
		      coeffs_dx->coeff_0_0 = t37;
		      coeffs_dx->coeff_p1_0 = t36+t40;
		      coeffs_dx->coeff_m1_p1 = t39+t41;
		      coeffs_dx->coeff_0_p1 = t37;
		      coeffs_dx->coeff_p1_p1 = t40+t42;

	    As illustrated, the code may use the macro RATIONAL
	    (defined later in this file) to represent rational-number
	    coefficients, and it may also declare temporary variables
	    as needed.

	    Note that this is the only included code which depends
	    on the actual interpolation scheme used; all the rest
	    just depends on the interpolation dimension and molecule
	    family and size.
  @endvar

  @version   $Header$
  @@*/

/******************************************************************************/

/*
 * ***** prototypes for functions local to this file *****
 */

#ifdef HAVE_OP_I
  static
    void compute_coeffs_I(FP_XYZ, struct COEFFS_STRUCT *coeffs_I);
#endif
#ifdef HAVE_OP_DX
  static
    void compute_coeffs_dx(FP_XYZ, struct COEFFS_STRUCT *coeffs_dx);
#endif
#ifdef HAVE_OP_DY
  static
    void compute_coeffs_dy(FP_XYZ, struct COEFFS_STRUCT *coeffs_dy);
#endif
#ifdef HAVE_OP_DZ
  static
    void compute_coeffs_dz(FP_XYZ, struct COEFFS_STRUCT *coeffs_dz);
#endif
#ifdef HAVE_OP_DXX
  static
    void compute_coeffs_dxx(FP_XYZ, struct COEFFS_STRUCT *coeffs_dxx);
#endif
#ifdef HAVE_OP_DXY
  static
    void compute_coeffs_dxy(FP_XYZ, struct COEFFS_STRUCT *coeffs_dxy);
#endif
#ifdef HAVE_OP_DXZ
  static
    void compute_coeffs_dxz(FP_XYZ, struct COEFFS_STRUCT *coeffs_dxz);
#endif
#ifdef HAVE_OP_DYY
  static
    void compute_coeffs_dyy(FP_XYZ, struct COEFFS_STRUCT *coeffs_dyy);
#endif
#ifdef HAVE_OP_DYZ
  static
    void compute_coeffs_dyz(FP_XYZ, struct COEFFS_STRUCT *coeffs_dyz);
#endif
#ifdef HAVE_OP_DZZ
  static
    void compute_coeffs_dzz(FP_XYZ, struct COEFFS_STRUCT *coeffs_dzz);
#endif

/******************************************************************************/

/*@@
  @routine	FUNCTION_NAME
  @date		23.Oct.2001
  @author	Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc
	This function does generalized interpolation of one or more
	2d arrays to arbitrary points.  For details, see the header
	comments for InterpLocalUniform() (in "InterpLocalUniform.c"
	in this same directory).

	This function's arguments are mostly all a subset of those of
	InterpLocalUniform() ; the difference is that this function
	takes all its arguments explicitly, whereas  InputLocalArrays()
	takes some of them indirectly via a key/value parameter table.
	Here we document only the "new" explicit arguments, and these
	only briefly; see the  InterpLocalUniform()  documentation
	and/or the thorn guide for this thorn for further details.

	If you change the arguments for this function, note that you
	must also change the prototype in "template.h".

  @var		log_fp
  @vdesc	This is NULL for no logging, or a valid stdio stream pointer
		to enable logging of the grid and interpolation coordinates.
  @vtype	FILE*
  @vio		in

  @var		error_info
  @vdesc	If error_info->per_point_status is non-NULL,
		   then we store per-point interpolator status for
			each point in  error_info->per_point_status[pt] .
		We set error_info->error_count to the number of
		points in error.
		We set the remaining elements of  error_info
		to a detailed description of the first point in error
		(if any).
  @vtype	struct error_info* error_info;
  @vio		out

  @var		molecule_structure_flags
  @vdesc	If this pointer is non-NULL, we store flags describing
		the interpolation molecule's structure in the pointed-to
		structure.
  @vtype	struct molecule_structure_flags* molecule_structure_flags;
  @vio		out

  @var		molecule_min_max_m_info
  @vdesc	If this pointer is non-NULL, we store the interpolation
		molecule's min/max m in the pointed-to structure.
  @vtype	struct molecule_min_max_m_info* molecule_min_max_m_info;
  @vio		out

  @var		molecule_positions
  @vdesc	If this pointer is non-NULL, we store the interpolation
		molecule's positions in the (caller-supplied) arrays
		pointed to by this pointer.
  @vtype	CCTK_INT* const molecule_positions[];
  @vio		out

  @var		Jacobian_info
  @vdesc	If this pointer is non-NULL, we store the Jacobian of
		the interpolation in the arrays (and in the manner)
		pointed to by this structure.
  @vtype	struct Jacobian_info* Jacobian_info;
  @vio		out

  @returntype   int
  @returndesc	This function's return results are a subset of those of
		InterpLocalUniform():
		0			success
		UTIL_ERROR_BAD_INPUT	one of the input arguments is invalid
		CCTK_ERROR_INTERP_POINT_OUTSIDE
					interpolation point is out of range
  @endreturndesc

  @@*/
int FUNCTION_NAME(/***** coordinate system *****/
		  const CCTK_REAL coord_origin[],
		  const CCTK_REAL coord_delta[],
		  /***** interpolation points *****/
		  int N_interp_points,
		  int interp_coords_type_code,
		  const void* const interp_coords[],
		  const CCTK_INT N_boundary_points_to_omit[],
		  const CCTK_REAL boundary_off_centering_tolerance[],
		  const CCTK_REAL boundary_extrapolation_tolerance[],
		  /***** input arrays *****/
		  int N_input_arrays,
		  const CCTK_INT input_array_offsets[],
		  const CCTK_INT input_array_strides[],
		  const CCTK_INT input_array_min_subscripts[],
		  const CCTK_INT input_array_max_subscripts[],
		  const CCTK_INT input_array_type_codes[],
		  const void* const input_arrays[],
		  /***** output arrays *****/
		  int N_output_arrays,
		  const CCTK_INT output_array_type_codes[],
		  void* const output_arrays[],
		  /***** operation info *****/
		  const CCTK_INT operand_indices[],
		  const CCTK_INT operation_codes[],
		  /***** debugging *****/
		  int debug, FILE* log_fp,
		  /***** other return results *****/
		  struct error_info* error_info,
		  struct molecule_structure_flags* molecule_structure_flags,
		  struct molecule_min_max_m_info* molecule_min_max_m_info,
		  CCTK_INT* const molecule_positions[],
		  struct Jacobian_info* Jacobian_info)
{
/*
 * ***** Naming conventions *****
 *
 * in, out = 0-origin indices each selecting an input/output array
 * pt = 0-origin index selecting an interpolation point
 * part = 0-origin index selecting real/imaginary part of a complex number
 */

/*
 * ***** Implementation notes: *****
 * 
 * The basic outline of this function is as follows:
 *
 *  store molecule structure flags
 *  if (querying molecule min/max m)
 *     then store molecule min/max m
 *  compute "which derivatives are wanted" flags
 *  precompute 1/dx factors
 *  set up input array strides
 *  set up Jacobian strides (or dummy values if no Jacobian info)
 *	for (int pt = 0 ; pt < N_interp_point ; ++pt)
 *	{
 *	struct DATA_STRUCT data;
 *	struct COEFFS_STRUCT coeffs_I, coeffs_dx, ..., coeffs_dzz;
 *	    for (int axis = 0 ; axis < N_dims ; ++axis)
 *	    {
 *	    switch  (interp_coords_type_codes)
 *		    {
 *	    case CCTK_VARIABLE_REAL:
 *		    load interp coords for this axis from
 *			 const CCTK_REAL *const interp_coords[]
 *		    break;
 *	    ...
 *		    }
 *	    }
 *	compute position of interpolation molecules
 *	   with respect to the grid into the XYZ variables
 *	   and store per-point status if requested
 *	if (this point is outside the grid)
 *	   then continue;
 *	if (querying molecule positions)
 *	   then store this molecule position
 *	compute molecule center 1-d position in input arrays
 *
 *	/# compute interpolation coefficients at this point #/
 *	if (want_I)
 *	   then compute_coeffs_I(XYZ, &coeffs_I);
 *	if (want_dx)
 *	   then compute_coeffs_dx(XYZ, &coeffs_dx);
 *	...
 *	if (want_dzz)
 *	   then compute_coeffs_dzz(XYZ, &coeffs_dx);
 *
 *	    for (int out = 0 ; out < N_output_arrays ; ++out)
 *	    {
 *	    const int in = operand_indices[out];
 *	    ***decode*** the output array datatype
 *	       to determine whether it's real or complex
 *	    const int N_output_parts = output array is complex ? 2 : 1;
 *		for (int part = 0 ; part < N_output_parts ; ++part)
 *		{
 *		if ( (input_arrays[in] != NULL)
 *		     && ( (input_arrays[in] != value at last load)
 *			  || (part != value at last load) ) )
 *		   then {
 *			save input_arrays[in] and part for
 *			   "previous value" test above
 *			***decode*** the input array datatype
 *			   to determine whether it's real or complex
 *			const int N_input_parts
 *			   = input array is complex ? 2 : 1;
 *			if (N_input_parts != N_output_parts)
 *			   then error(...)
 *			switch  (input_array_type_codes[in])
 *				{
 *			case CCTK_VARIABLE_REAL:
 *				LOAD_DATA_REAL(input_array_ptr,
 *					       STRIDE_IJK,
 *					       &data);
 *				break;
 *			    ...
 *			case CCTK_VARIABLE_COMPLEX:
 *				LOAD_DATA_COMPLEX(input_array_ptr,
 *						  STRIDE_IJK, part,
 *						  &data);
 *				break;
 *			    ...
 *				}
 *			}
 *
 *		if (output_arrays[out] != NULL)
 *		   then {
 *			/# interpolate at this point #/
 *			fp result;
 *			switch  (operation_codes[out])
 *				{
 *			case 0:
 *				result = EVALUATE_MOLECULE(&coeffs_I,
 *							   &data);
 *				break;
 *			case 1:
 *				result = inverse_dx
 *					 * EVALUATE_MOLECULE(&coeffs_dx,
 *							     &data);
 *				break;
 *			...
 *			case 33:
 *				result = inverse_dz * inverse_dz
 *					 * EVALUATE_MOLECULE(&coeffs_dzz,
 *							     &data);
 *				break;
 *				}
 *
 *			/# store result in output array #/
 *			switch  (output_array_type_codes[out])
 *				{
 *			case CCTK_VARIABLE_REAL:
 *				store result in
 *				       CCTK_REAL *const output_arrays[]
 *				break;
 *			...
 *			case CCTK_VARIABLE_COMPLEX:
 *				store result in
 *				       CCTK_COMPLEX *const output_arrays[][part]
 *				break;
 *			...
 *				}
 *			}
 *		if (querying Jacobian && (Jacobian_pointer[out] != NULL))
 *		   then {
 *			/# store Jacobian coefficients at this point #/
 *			CCTK_REAL *const Jacobian_ptr
 *				= Jacobian_info->Jacobian_pointer[out]
 *				  + Jacobian_info->Jacobian_offset[out];
 *			switch	(operation_code)
 *				{
 *			case 0:
 *				STORE_COEFFS(1.0, &coeffs_I,
 *					     Jacobian_ptr, JACOBIAN_MIJK_STRIDE,
 *					     pt, part,
 *				break;
 *			case 1:
 *				STORE_COEFFS(inverse_dx, &coeffs_dx,
 *					     Jacobian_ptr, JACOBIAN_MIJK_STRIDE,
 *					     pt, part,
 *				break;
 *			...
 *			case 33:
 *				STORE_COEFFS(inverse_dz*inverse_dz, &coeffs_dzz,
 *					     Jacobian_ptr, JACOBIAN_MIJK_STRIDE,
 *					     pt, part,
 *				break;
 *				}
 *			}
 *		}
 *	    }
 *	}
 *
 * For complex datatypes we pointer-alias the N-dimensional input array
 * to a 2-dimensional array where the 1st axis is the 1-D subscript
 * corresponding to the N input axes, and the 2nd axes has 2 elements
 * subscripted by [part] for the (real,imaginary) components of the
 * complex values.
 *
 * We do all floating-point computations in type "fp" (typically a typedef
 * for CCTK_REAL), so arrays of higher precision than this will incur extra
 * rounding errors.  In practice these should be negligible compared to the
 * interpolation errors.
 */


/* layout of axes in N_dims-element arrays */
#define X_AXIS	0
#define Y_AXIS	1
#define Z_AXIS	2
#if (N_DIMS > 3)
  #error "N_DIMS may not be > 3!"
#endif

/* layout of axes and min/max ends in out_of_range_tolerance[] array */
#define X_AXIS_MIN	0
#define X_AXIS_MAX	1
#define Y_AXIS_MIN	2
#define Y_AXIS_MAX	3
#define Z_AXIS_MIN	4
#define Z_AXIS_MAX	5
#if (N_DIMS > 3)
  #error "N_DIMS may not be > 3!"
#endif

/* basic sanity check on molecule size */
#define MOLECULE_M_COUNT (MOLECULE_MAX_M - MOLECULE_MIN_M + 1)
#if (MOLECULE_SIZE != MOLECULE_M_COUNT)
  #error "MOLECULE_SIZE inconsistent with MOLECULE_{MIN,MAX}_M!"
#endif


/* macros used by machine-generated interpolation coefficient expressions */
/*
 * FIXME: right now this is used as (eg) RATIONAL(1.0,2.0);
 *	  it might be cleaner if it were RATIONAL(1,2) with the
 *	  preprocessor ## operator used to glue on the .0
 *	  (I _think_ this is portable, but is it really?)
 */
#define RATIONAL(num,den)	(num/den)


/*
 * store molecule structure flags, molecule min/max m (if requested)
 */
if (molecule_structure_flags != NULL)
   then {
	molecule_structure_flags->MSS_is_fn_of_interp_coords = 0;
	molecule_structure_flags->MSS_is_fn_of_which_operation = 0;
	molecule_structure_flags->MSS_is_fn_of_input_array_values = 0;
	molecule_structure_flags->Jacobian_is_fn_of_input_array_values = 0;
	}
if (molecule_min_max_m_info != NULL)
   then {
	int axis;
		for (axis = 0 ; axis < N_DIMS ; ++axis)
		{
		molecule_min_max_m_info->molecule_min_m[axis] = MOLECULE_MIN_M;
		molecule_min_max_m_info->molecule_max_m[axis] = MOLECULE_MAX_M;
		}
	}


/*
 * compute "which derivatives are wanted" flags
 */
  {
#ifdef HAVE_OP_I
  bool want_I = false;
#endif
#ifdef HAVE_OP_DX
  bool want_dx = false;
#endif
#ifdef HAVE_OP_DY
  bool want_dy = false;
#endif
#ifdef HAVE_OP_DZ
  bool want_dz = false;
#endif
#ifdef HAVE_OP_DXX
  bool want_dxx = false;
#endif
#ifdef HAVE_OP_DXY
  bool want_dxy = false;
#endif
#ifdef HAVE_OP_DXZ
  bool want_dxz = false;
#endif
#ifdef HAVE_OP_DYY
  bool want_dyy = false;
#endif
#ifdef HAVE_OP_DYZ
  bool want_dyz = false;
#endif
#ifdef HAVE_OP_DZZ
  bool want_dzz = false;
#endif

  {
int out;
	for (out = 0 ; out < N_output_arrays ; ++out)
	{
	switch	(operation_codes[out])
		{
	#ifdef HAVE_OP_I
	  case 0:	want_I = true;		break;
	#endif
	#ifdef HAVE_OP_DX
	  case 1:	want_dx = true;		break;
	#endif
	#ifdef HAVE_OP_DY
	  case 2:	want_dy = true;		break;
	#endif
	#ifdef HAVE_OP_DZ
	  case 3:	want_dz = true;		break;
	#endif
	#ifdef HAVE_OP_DXX
	  case 11:	want_dxx = true;	break;
	#endif
	#ifdef HAVE_OP_DXY
	  case 12:
	  case 21:	want_dxy = true;	break;
	#endif
	#ifdef HAVE_OP_DXZ
	  case 13:
	  case 31:	want_dxz = true;	break;
	#endif
	#ifdef HAVE_OP_DYY
	  case 22:	want_dyy = true;	break;
	#endif
	#ifdef HAVE_OP_DYZ
	  case 23:
	  case 32:	want_dyz = true;	break;
	#endif
	#ifdef HAVE_OP_DZZ
	  case 33:	want_dzz = true;	break;
	#endif
	default:
	  CCTK_CVWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): generalized interpolation operation_code %d\n"
"                              is not supported!\n"
"                              (0-origin) output #out=%d"
		     ,
		     (int) operation_codes[out],
		     out);					/*NOTREACHED*/
		return UTIL_ERROR_BAD_INPUT;		/*** ERROR RETURN ***/
		}
	}
  }


/*
 * save origin/delta variables, precompute 1/delta factors
 * ... in theory we could compute only those factors we're going to use,
 *     but it's not worth the trouble, so we just compute them all
 */
  {
#if N_DIMS >= 1
  #if    defined(HAVE_OP_DX) \
      || defined(HAVE_OP_DXY) || defined(HAVE_OP_DXZ) \
      || defined(HAVE_OP_DXX)
    const fp inverse_dx = 1.0 / coord_delta[X_AXIS];
  #endif
#endif
#if N_DIMS >= 2
  #if    defined(HAVE_OP_DY) \
      || defined(HAVE_OP_DXY) || defined(HAVE_OP_DYZ) \
      || defined(HAVE_OP_DYY)
    const fp inverse_dy = 1.0 / coord_delta[Y_AXIS];
  #endif
#endif
#if N_DIMS >= 3
  #if    defined(HAVE_OP_DZ) \
      || defined(HAVE_OP_DXZ) || defined(HAVE_OP_DYZ) \
      || defined(HAVE_OP_DZZ)
    const fp inverse_dz = 1.0 / coord_delta[Z_AXIS];
  #endif
#endif
#if (N_DIMS > 3)
  #error "N_DIMS may not be > 3!"
#endif


/*
 * set up input array strides
 */
#if (N_DIMS >= 1)
  const int stride_i = input_array_strides[X_AXIS];
#endif
#if (N_DIMS >= 2)
  const int stride_j = input_array_strides[Y_AXIS];
#endif
#if (N_DIMS >= 3)
  const int stride_k = input_array_strides[Z_AXIS];
#endif
#if (N_DIMS > 3)
  #error "N_DIMS may not be > 3!"
#endif


/*
 * set up Jacobian m strides, or dummy values if no Jacobian info
 */
#if N_DIMS >= 1
  const int Jacobian_mi_stride = (Jacobian_info == NULL)
				 ? 0
				 : Jacobian_info->Jacobian_m_strides[X_AXIS];
#endif
#if (N_DIMS >= 2)
  const int Jacobian_mj_stride = (Jacobian_info == NULL)
				 ? 0
				 : Jacobian_info->Jacobian_m_strides[Y_AXIS];
#endif
#if (N_DIMS >= 3)
  const int Jacobian_mk_stride = (Jacobian_info == NULL)
				 ? 0
				 : Jacobian_info->Jacobian_m_strides[Z_AXIS];
#endif
#if (N_DIMS > 3)
  #error "N_DIMS may not be > 3!"
#endif


/*
 * if requested, write the interpolation grid to the log file
 * (we'll write (some of) the individual interpolation coordinates
 *  as we interpolate them (below))
 */
if (log_fp != NULL)
   then {
	fprintf(log_fp, "%d\t%d", (int) N_DIMS, N_interp_points);

	  {
	fp grid_min_xyz[MAX_N_DIMS], grid_max_xyz[MAX_N_DIMS];
	int axis;
	    for (axis = 0 ; axis < N_DIMS ; ++axis)
	    {
	    grid_min_xyz[axis]
		= coord_origin[axis]
		  + input_array_min_subscripts[axis]*coord_delta[axis];
	    grid_max_xyz[axis]
		= coord_origin[axis]
		  + input_array_max_subscripts[axis]*coord_delta[axis];
	    }

      #if   (N_DIMS == 1)
	fprintf(log_fp, "\t%g\t%g\t%g\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0",
		(double) grid_min_xyz[X_AXIS],
		(double) coord_delta [X_AXIS],
		(double) grid_max_xyz[X_AXIS]);
      #elif (N_DIMS == 2)
	fprintf(log_fp, "\t%g\t%g\t%g\t%g\t%g\t%g\t0.0\t0.0\t0.0",
		(double) grid_min_xyz[X_AXIS],
		(double) coord_delta [X_AXIS],
		(double) grid_max_xyz[X_AXIS],
		(double) grid_min_xyz[Y_AXIS],
		(double) coord_delta [Y_AXIS],
		(double) grid_max_xyz[Y_AXIS]);
      #elif (N_DIMS == 3)
	fprintf(log_fp, "\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g",
		(double) grid_min_xyz[X_AXIS],
		(double) coord_delta [X_AXIS],
		(double) grid_max_xyz[X_AXIS],
		(double) grid_min_xyz[Y_AXIS],
		(double) coord_delta [Y_AXIS],
		(double) grid_max_xyz[Y_AXIS],
		(double) grid_min_xyz[Z_AXIS],
		(double) coord_delta [Z_AXIS],
		(double) grid_max_xyz[Z_AXIS]);
      #else
	#error "N_DIMS may not be > 3!"
      #endif
	  }
	}


/*
 * interpolate at each point
 */
  {
int pt;
int return_status = 0;
error_info->error_count = 0;

if (debug > 0)
   then {
	printf("AEILocalInterp:: in interpolator fn: N_DIMS=%d N_interp_points=%d\n", (int) N_DIMS, N_interp_points);
	fflush(stdout);
	}

    for (pt = 0 ; pt < N_interp_points ; ++pt)
    {
    /* this struct holds a molecule-sized piece of a single */
    /* real input array, or of a real/complex component of */
    /* a complex input array */
    struct DATA_STRUCT data;

    /* each of these structs holds a molecule-sized set of */
    /* interpolation coefficients -- recall we are representing */
    /* the interpolant as a linear combination of the data values */
    #ifdef HAVE_OP_I
      struct COEFFS_STRUCT coeffs_I;
    #endif
    #ifdef HAVE_OP_DX
      struct COEFFS_STRUCT coeffs_dx;
    #endif
    #ifdef HAVE_OP_DY
      struct COEFFS_STRUCT coeffs_dy;
    #endif
    #ifdef HAVE_OP_DZ
      struct COEFFS_STRUCT coeffs_dz;
    #endif
    #ifdef HAVE_OP_DXX
      struct COEFFS_STRUCT coeffs_dxx;
    #endif
    #ifdef HAVE_OP_DXY
      struct COEFFS_STRUCT coeffs_dxy;
    #endif
    #ifdef HAVE_OP_DXZ
      struct COEFFS_STRUCT coeffs_dxz;
    #endif
    #ifdef HAVE_OP_DYY
      struct COEFFS_STRUCT coeffs_dyy;
    #endif
    #ifdef HAVE_OP_DYZ
      struct COEFFS_STRUCT coeffs_dyz;
    #endif
    #ifdef HAVE_OP_DZZ
      struct COEFFS_STRUCT coeffs_dzz;
    #endif


    /*
     * load the interpolation point coordinates
     * from the interp_coords[] arrays
     * FIXME: Maybe it would be better (faster) to do this
     *	  with N_DIMS open-coded calls on a function?
     *	  But then we'd have to have a sentinal value
     *	  return for the unknown-type-code error case.
     *	  Yuk! :( :(
     */
    fp interp_coords_fp[N_DIMS];
      {
    int axis;
	for (axis = 0 ; axis < N_DIMS ; ++axis)
	{
	/* pointer to array of interp coords for this axis */
	const void *const interp_coords_ptr = interp_coords[axis];

	switch	(interp_coords_type_code)
		{
	case CCTK_VARIABLE_REAL:
		  {
		const CCTK_REAL *const interp_coords_ptr_real
			= (const CCTK_REAL *) interp_coords_ptr;
		interp_coords_fp[axis] = interp_coords_ptr_real[pt];
		break;
		  }

	#ifdef HAVE_CCTK_REAL4
	case CCTK_VARIABLE_REAL4:
		  {
		const CCTK_REAL4 *const interp_coords_ptr_real4
			= (const CCTK_REAL4 *) interp_coords_ptr;
		interp_coords_fp[axis] = interp_coords_ptr_real4[pt];
		break;
		  }
	#endif

	#ifdef HAVE_CCTK_REAL8
	case CCTK_VARIABLE_REAL8:
		  {
		const CCTK_REAL8 *const interp_coords_ptr_real8
			= (const CCTK_REAL8 *) interp_coords_ptr;
		interp_coords_fp[axis] = interp_coords_ptr_real8[pt];
		break;
		  }
	#endif

	#ifdef HAVE_CCTK_REAL16
	case CCTK_VARIABLE_REAL16:
		  {
		/* FIXME: maybe we should warn (once per cactus run) */
		/*        that we may be doing arithmetic in lower */
		/*        precision than the interp coords? */
		const CCTK_REAL16 *const interp_coords_ptr_real16
			= (const CCTK_REAL16 *) interp_coords_ptr;
		interp_coords_fp[axis]
			= interp_coords_ptr_real16[pt];
		break;
		  }
	#endif

	default:
		CCTK_CVWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): interpolation coordinates datatype %d\n"
"                              is not supported!"
			   ,
			   interp_coords_type_code);
							/*NOTREACHED*/
		return UTIL_ERROR_BAD_INPUT;	/*** ERROR RETURN ***/
		/* end of switch (interp_coords_type_code) */
		}

	/* end of for (axis = ...) loop to get interp coords */
	}
      }

    if (debug >= 6)
       then {
      #if   (N_DIMS == 1)
	    printf("pt=%d interp coord %.16g\n",
		   pt, (double) interp_coords_fp[X_AXIS]);
      #elif (N_DIMS == 2)
	    printf("pt=%d interp coords %.16g %.16g\n",
		   pt, (double) interp_coords_fp[X_AXIS],
		       (double) interp_coords_fp[Y_AXIS]);
      #elif (N_DIMS == 3)
	    printf("pt=%d interp coords %.16g %.16g %.16g\n",
		   pt, (double) interp_coords_fp[X_AXIS],
		       (double) interp_coords_fp[Y_AXIS],
		       (double) interp_coords_fp[Z_AXIS]);
      #else
	    #error "N_DIMS may not be > 3!"
      #endif
	    fflush(stdout);
	    }

    /*
     * if requested, write this point's interpolation coords to the log file
     */
    if ( (log_fp != NULL)
	 && ((pt < 3) || (pt == N_interp_points-1)) )
       then {
      #if   (N_DIMS == 1)
	    fprintf(log_fp, "\t%g\t0.0\t0.0",
		    (double) interp_coords_fp[X_AXIS]);
      #elif (N_DIMS == 2)
	    fprintf(log_fp, "\t%g\t%g\t0.0",
		    (double) interp_coords_fp[X_AXIS],
		    (double) interp_coords_fp[Y_AXIS]);
      #elif (N_DIMS == 3)
	    fprintf(log_fp, "\t%g\t%g\t%g",
		    (double) interp_coords_fp[X_AXIS],
		    (double) interp_coords_fp[Y_AXIS],
		    (double) interp_coords_fp[Z_AXIS]);
      #else
	#error "N_DIMS may not be > 3!"
      #endif
	    }

    /*
     * compute position of interpolation molecules with respect to
     * the grid, i.e. compute (x,y,z) coordinates of interpolation
     * point relative to molecule center, in units of the grid spacing
     *
     * n.b. we need the final answers in variables with the magic
     *      names (x,y,z) (machine-generated code uses these names),
     *      but we use temp variables as intermediates for these and
     *      for  center_[ijk] for (likely) better performance:
     *      the temp variables have their addresses taken and so
     *      may not be register-allocated, whereas the final variables
     *	    are "clean" and thus more likely to be register-allocated
     */
      {
    int this_point_status = 0;
    int center_ijk_temp[MAX_N_DIMS];
    fp  xyz_temp[MAX_N_DIMS];
    int axis;
	for (axis = 0 ; axis < N_DIMS ; ++axis)
	{
	const int ibndry_min = 2*axis;
	const int ibndry_max = 2*axis + 1;
	if (debug >= 8)
	   then {
		printf("axis=%d   ibndry_{min,max}=%d,%d   MOLECULE_SIZE=%d\n",
		       axis, ibndry_min, ibndry_max, MOLECULE_SIZE);
		printf("   coord_origin[%d]=%g coord_delta[%d]=%g\n",
		       axis, (double)coord_origin[axis],
		       axis, (double)coord_delta[axis]);
		printf("   input_array_[min,max]_subscripts[%d]=[%d,%d]\n",
		       axis,
		       (int)input_array_min_subscripts[axis],
		       (int)input_array_max_subscripts[axis]);
		printf("   N_boundary_points_to_omit[ibndry_[min,max]]=[%d,%d]\n",
		       (int)N_boundary_points_to_omit[ibndry_min],
		       (int)N_boundary_points_to_omit[ibndry_max]);
		}
	  {
	const int mp_status = AEILocalInterp_molecule_posn
				 (coord_origin[axis], coord_delta[axis],
				  input_array_min_subscripts[axis]
				    + N_boundary_points_to_omit[ibndry_min],
				  input_array_max_subscripts[axis]
				    - N_boundary_points_to_omit[ibndry_max],
				  MOLECULE_SIZE,
				  boundary_off_centering_tolerance[ibndry_min],
				  boundary_off_centering_tolerance[ibndry_max],
				  boundary_extrapolation_tolerance[ibndry_min],
				  boundary_extrapolation_tolerance[ibndry_max],
				  interp_coords_fp[axis],
				  debug,
				  &center_ijk_temp[axis], &xyz_temp[axis]);
	if (debug >= 8)
	   then printf("==> got mp_status=%d\n", mp_status);

	/*
	 * unfortunately, the status codes from AEILocalInterp_molecule_posn()
	 * are different from the ones we need ==> we have to decode them
	 * to get the status for this axis
	 */
	  {

	int this_axis_status = -1;
	switch	(mp_status)
		{
	case 0:
		this_axis_status = 0;
		break;
	case MOLECULE_POSN_ERROR_GRID_TINY:
		this_axis_status = CCTK_ERROR_INTERP_GRID_TOO_SMALL;
		break;
	case MOLECULE_POSN_ERROR_X_LT_MIN:
	case MOLECULE_POSN_ERROR_X_GT_MAX:
		this_axis_status = CCTK_ERROR_INTERP_POINT_OUTSIDE;
		break;
	case MOLECULE_POSN_ERROR_NAN:
		this_axis_status = CCTK_ERROR_INTERP_COORD_NAN;
		break;
	case MOLECULE_POSN_ERROR_DX_ZERO:
		this_axis_status = CCTK_ERROR_INTERP_DELTA_X_ZERO;
		break;
	default:
		CCTK_CVWarn(BUG_MSG_SEVERITY_LEVEL,
			   __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): internal error!\n"
"        unexpected result %d from AEILocalInterp_molecule_posn()\n"
"        at pt=%d axis=%d!\n"
			   ,
			   mp_status, pt, axis);		/*NOTREACHED*/
		}

	if (this_axis_status < this_point_status)
	   then this_point_status = this_axis_status;

	/* end of for (axis = ...) loop to locate interp points in the grid */
	  }
	  }
	}

    if (error_info->per_point_status != NULL)
       then error_info->per_point_status[pt] = this_point_status;

    if (this_point_status < 0)
       then {
	    ++error_info->error_count;

	    if (error_info->print_warning_msg)
	       then {
		    fp grid_min_xyz[MAX_N_DIMS], grid_max_xyz[MAX_N_DIMS];
			for (axis = 0 ; axis < N_DIMS ; ++axis)
			{
			grid_min_xyz[axis]
			 = coord_origin[axis]
			   + input_array_min_subscripts[axis]*coord_delta[axis];
			grid_max_xyz[axis]
			 = coord_origin[axis]
			   + input_array_max_subscripts[axis]*coord_delta[axis];
			}
		    CCTK_CVWarn(WARNING_MSG_SEVERITY_LEVEL,
			       __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform():\n"
"        interpolation point is either outside the grid,\n"
"        or inside but too close to the grid boundary!\n"
"        (this may be caused by a global interpolation with\n"
"         driver::ghost_size too small)\n"
"        0-origin interpolation point number pt=%d of N_interp_points=%d\n"
#if   (N_DIMS == 1)
"        interpolation point x=%g\n"
"        grid x_min(delta_x)x_max = %g(%g)%g\n"
			       ,
			       pt, N_interp_points,
	   (double) interp_coords_fp[X_AXIS],
	   (double) grid_min_xyz[X_AXIS], (double) coord_delta [X_AXIS],
					  (double) grid_max_xyz[X_AXIS]);
#elif (N_DIMS == 2)
"        interpolation point (x,y)=(%g,%g)\n"
"        grid x_min(delta_x)x_max = %g(%g)%g\n"
"        grid y_min(delta_y)y_max = %g(%g)%g\n"
			       ,
			       pt, N_interp_points,
	   (double) interp_coords_fp[X_AXIS], (double) interp_coords_fp[Y_AXIS],
	   (double) grid_min_xyz[X_AXIS], (double) coord_delta [X_AXIS],
					  (double) grid_max_xyz[X_AXIS],
	   (double) grid_min_xyz[Y_AXIS], (double) coord_delta [Y_AXIS],
					  (double) grid_max_xyz[Y_AXIS]);
#elif (N_DIMS == 3)
"        interpolation point (x,y,z)=(%g,%g,%g)\n"
"        grid x_min(delta_x)x_max = %g(%g)%g\n"
"        grid y_min(delta_y)y_max = %g(%g)%g\n"
"        grid z_min(delta_z)z_max = %g(%g)%g\n"
			       ,
			       pt, N_interp_points,
	   (double) interp_coords_fp[X_AXIS], (double) interp_coords_fp[Y_AXIS],
					      (double) interp_coords_fp[Z_AXIS],
	   (double) grid_min_xyz[X_AXIS], (double) coord_delta [X_AXIS],
					  (double) grid_max_xyz[X_AXIS],
	   (double) grid_min_xyz[Y_AXIS], (double) coord_delta [Y_AXIS],
					  (double) grid_max_xyz[Y_AXIS],
	   (double) grid_min_xyz[Z_AXIS], (double) coord_delta [Z_AXIS],
					  (double) grid_max_xyz[Z_AXIS]);
#else
  #error "N_DIMS may not be > 3!"
#endif
		    }

	    if (this_point_status < return_status)
	       then return_status = this_point_status;

	    if (debug >= 6)
	       then {
		    printf("AEILocalInterp:: pt=%d has this_point_status=%d ==> skipping it\n", pt, this_point_status);
		    fflush(stdout);
		    }

	    /* this point is in error (eg outside grid) */
	    /* ==> don't try to interpolate at this point! */
	    continue;
	    }

      {
  #if (N_DIMS >= 1)
    const int center_i = center_ijk_temp[X_AXIS];
    const fp  x        = xyz_temp       [X_AXIS];
  #endif
  #if (N_DIMS >= 2)
    const int center_j = center_ijk_temp[Y_AXIS];
    const fp  y        = xyz_temp       [Y_AXIS];
  #endif
  #if (N_DIMS >= 3)
    const int center_k = center_ijk_temp[Z_AXIS];
    const fp  z        = xyz_temp       [Z_AXIS];
  #endif

    if (debug >= 6)
       then {
	  #if   (N_DIMS == 1)
	    printf("interp center_i = %d\n", center_i);
	    printf("interp grid-relative x = %.16g\n", (double) x);
	  #elif (N_DIMS == 2)
	    printf("interp center_ij = %d %d\n", center_i, center_j);
	    printf("interp grid-relative xy = %.16g %.16g\n",
		   (double) x, (double) y);
	  #elif (N_DIMS == 3)
	    printf("interp center_ijk = %d %d %d\n",
		   center_i, center_j, center_k);
	    printf("interp grid-relative xyz = %.16g %.16g %.16g\n",
		   (double) x, (double) y, (double) z);
	  #else
	    #error "N_DIMS may not be > 3!"
	  #endif
	    fflush(stdout);
	    }

    /*
     * compute 1-d position of molecule center in input arrays
     */
      {
    #if   (N_DIMS == 1)
      const int molecule_center_posn =   stride_i*center_i;
    #elif (N_DIMS == 2)
      const int molecule_center_posn =   stride_i*center_i
				       + stride_j*center_j;
    #elif (N_DIMS == 3)
      const int molecule_center_posn =   stride_i*center_i
				       + stride_j*center_j
				       + stride_k*center_k;
    #else
      #error "N_DIMS may not be > 3!"
    #endif


    /*
     * molecule position queries
     */
    if (molecule_positions != NULL)
       then {
	    #if (N_DIMS >= 1)
	      molecule_positions[X_AXIS][pt] = center_i;
	    #endif
	    #if (N_DIMS >= 2)
	      molecule_positions[Y_AXIS][pt] = center_j;
	    #endif
	    #if (N_DIMS >= 3)
	      molecule_positions[Z_AXIS][pt] = center_k;
	    #endif
	    }


    /*
     * compute the coefficients at this point for whichever
     * operation_codes[] values (derivatatives) are wanted
     * (these are polynomials in the variables (x,y,z))
     */
    #ifdef HAVE_OP_I
      if (want_I)
	 then compute_coeffs_I(XYZ, &coeffs_I);
    #endif
    #ifdef HAVE_OP_DX
      if (want_dx)
	 then compute_coeffs_dx(XYZ, &coeffs_dx);
    #endif
    #ifdef HAVE_OP_DY
      if (want_dy)
	 then compute_coeffs_dy(XYZ, &coeffs_dy);
    #endif
    #ifdef HAVE_OP_DZ
      if (want_dz)
	 then compute_coeffs_dz(XYZ, &coeffs_dz);
    #endif
    #ifdef HAVE_OP_DXX
      if (want_dxx)
	 then compute_coeffs_dxx(XYZ, &coeffs_dxx);
    #endif
    #ifdef HAVE_OP_DXY
      if (want_dxy)
	 then compute_coeffs_dxy(XYZ, &coeffs_dxy);
    #endif
    #ifdef HAVE_OP_DXZ
      if (want_dxz)
	 then compute_coeffs_dxz(XYZ, &coeffs_dxz);
    #endif
    #ifdef HAVE_OP_DYY
      if (want_dyy)
	 then compute_coeffs_dyy(XYZ, &coeffs_dyy);
    #endif
    #ifdef HAVE_OP_DYZ
      if (want_dyz)
	 then compute_coeffs_dyz(XYZ, &coeffs_dyz);
    #endif
    #ifdef HAVE_OP_DZZ
      if (want_dzz)
	 then compute_coeffs_dzz(XYZ, &coeffs_dzz);
    #endif


	/*
	 * compute each output array at this point
	 */
	  {
	int out;

	/*
	 * next 2 initializers must be invalid values to make sure we
	 * execute the ***load*** the first time in the test at the
	 * top of the  part  loop below
	 */
	const void* input_array_ptr__last_load = NULL;
	int part__last_load = -1;

	    for (out = 0 ; out < N_output_arrays ; ++out)
	    {
	    const int in = operand_indices[out];
	    const void* const input_array_ptr = input_arrays[in];

	    /*
 	     * ***decode*** the output array datatype
	     * to determine whether it's real or complex,
	     */
	    const int N_output_parts
	       = AEILocalInterp_decode_N_parts(output_array_type_codes[out]);
	    if (! ((N_output_parts == 1) || (N_output_parts == 2)))
	       then {
CCTK_CVWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform():\n"
"        output array doesn't seem to be a real or complex number,\n"
"        or more precisely, output array has number of \"real parts\"\n"
"        (1=real, 2=complex) which isn't 1 or 2!\n"
"        0-origin output #out=%d\n"
"        datatype code=%d\n"
"           (datatype codes are defined by the Cactus flesh,\n"
"            see src/include/cctk_Constants.h)\n"
"        ==> N_parts=%d"
	   ,
	   out, (int) output_array_type_codes[out], N_output_parts);
								/*NOTREACHED*/
		    return UTIL_ERROR_BAD_INPUT;	/*** ERROR RETURN ***/
		    }

	      {
	    int part;
		for (part = 0 ; part < N_output_parts ; ++part)
		{
		if ( (input_array_ptr != NULL)
		     &&
		     ((input_array_ptr != input_array_ptr__last_load)
		      || (part != part__last_load)) )
		   then {
			/* remember when we did the following load */
			input_array_ptr__last_load = input_array_ptr;
			part__last_load = part;

			/*
			 * ***decode*** the input array datatype
			 * to determine whether it's real or complex,
			 */
			  {
			const int N_input_parts
			   = AEILocalInterp_decode_N_parts(input_array_type_codes[in]);
			if (N_input_parts != N_output_parts)
			   then {
CCTK_CVWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform():\n"
"        data types are incompatible between input and output arrays, or\n"
"        more precisely, number of \"real parts\" (1=real, 2=complex) differ!\n"
"        0-origin input  #in =%d datatype code=%d N_parts=%d\n"
"        0-origin output #out=%d datatype code=%d N_parts=%d\n"
"           (datatype codes are defined by the Cactus flesh,\n"
"            see src/include/cctk_Constants.h)"
	   ,
	   in, (int) input_array_type_codes[in], N_input_parts,
	   out, (int) output_array_type_codes[out], N_output_parts);
								/*NOTREACHED*/
				return UTIL_ERROR_BAD_INPUT;
							/*** ERROR RETURN ***/
				}

			/*
			 * load the molecule-sized piece of
			 *  input_arrays[in][part]  at this molecule
			 * position, into the data struct
			 */
			  {
			const int input_posn = molecule_center_posn
					       + input_array_offsets[in];
			switch	(input_array_type_codes[in])
				{
case CCTK_VARIABLE_REAL:
	  {
	const CCTK_REAL *const input_array_ptr_real
		= ((const CCTK_REAL *) input_array_ptr) + input_posn;
	LOAD_DATA_REAL(input_array_ptr_real,
		       STRIDE_IJK,
		       &data);
      #if (N_DIMS == 2) && (MOLECULE_SIZE == 4)
	if (debug >= 10)
	   then {
		/* we assume a cubical molecule :( */
		printf("AEILocalInterp:: loaded data is:\n");
		printf("   data_m1_m1 = %g\n", (double) data.data_m1_m1);
		printf("   data_0_m1  = %g\n", (double) data.data_0_m1);
		printf("   data_p1_m1 = %g\n", (double) data.data_p1_m1);
		printf("   data_p2_m1 = %g\n", (double) data.data_p2_m1);
		printf("   data_m1_0  = %g\n", (double) data.data_m1_0);
		printf("   data_0_0   = %g\n", (double) data.data_0_0);
		printf("   data_p1_0  = %g\n", (double) data.data_p1_0);
		printf("   data_p2_0  = %g\n", (double) data.data_p2_0);
		printf("   data_m1_p1 = %g\n", (double) data.data_m1_p1);
		printf("   data_0_p1  = %g\n", (double) data.data_0_p1);
		printf("   data_p1_p1 = %g\n", (double) data.data_p1_p1);
		printf("   data_p2_p1 = %g\n", (double) data.data_p2_p1);
		printf("   data_m1_p2 = %g\n", (double) data.data_m1_p2);
		printf("   data_0_p2  = %g\n", (double) data.data_0_p2);
		printf("   data_p1_p2 = %g\n", (double) data.data_p1_p2);
		printf("   data_p2_p2 = %g\n", (double) data.data_p2_p2);
		fflush(stdout);
		}
      #endif
	break;
	  }

#ifdef HAVE_CCTK_REAL4
case CCTK_VARIABLE_REAL4:
	  {
	const CCTK_REAL4 *const input_array_ptr_real4
		= ((const CCTK_REAL4 *) input_array_ptr) + input_posn;
	LOAD_DATA_REAL4(input_array_ptr_real4,
			STRIDE_IJK,
			&data);
	break;
	  }
#endif

#ifdef HAVE_CCTK_REAL8
case CCTK_VARIABLE_REAL8:
	  {
	const CCTK_REAL8 *const input_array_ptr_real8
		= ((const CCTK_REAL8 *) input_array_ptr) + input_posn;
	LOAD_DATA_REAL8(input_array_ptr_real8,
			STRIDE_IJK,
			&data);
	break;
	  }
#endif

#ifdef HAVE_CCTK_REAL16
case CCTK_VARIABLE_REAL16:
	  {
	/* FIXME: maybe we should warn (once per cactus run) that we may be */
	/*        doing arithmetic in lower precision than the data type? */
	const CCTK_REAL16 *const input_array_ptr_real16
		= ((const CCTK_REAL16 *) input_array_ptr) + input_posn;
	LOAD_DATA_REAL16(input_array_ptr_real16,
			 STRIDE_IJK,
			 &data);
	break;
	  }
#endif

case CCTK_VARIABLE_COMPLEX:
	  {
	const CCTK_REAL (*const input_array_ptr_complex)[COMPLEX_N_PARTS]
		= ((const CCTK_REAL (*)[COMPLEX_N_PARTS]) input_array_ptr)
		  + input_posn;
	LOAD_DATA_COMPLEX(input_array_ptr_complex,
			  STRIDE_IJK, part,
			  &data);
	break;
	  }

#ifdef HAVE_CCTK_COMPLEX8
case CCTK_VARIABLE_COMPLEX8:
	  {
	const CCTK_REAL4 (*const input_array_ptr_complex8)[COMPLEX_N_PARTS]
		= ((const CCTK_REAL4 (*)[COMPLEX_N_PARTS]) input_array_ptr)
		  + input_posn;
	LOAD_DATA_COMPLEX8(input_array_ptr_complex8,
			   STRIDE_IJK, part,
			   &data);
	break;
	  }
#endif

#ifdef HAVE_CCTK_COMPLEX16
case CCTK_VARIABLE_COMPLEX16:
	  {
	const CCTK_REAL8 (*const input_array_ptr_complex16)[COMPLEX_N_PARTS]
		= ((const CCTK_REAL8 (*)[COMPLEX_N_PARTS]) input_array_ptr)
		  + input_posn;
	LOAD_DATA_COMPLEX16(input_array_ptr_complex16,
			    STRIDE_IJK, part,
			    &data);
	break;
	  }
#endif


#ifdef HAVE_CCTK_COMPLEX32
case CCTK_VARIABLE_COMPLEX32:
	  {
	/* FIXME: maybe we should warn (once per cactus run) that we may be */
	/*        doing arithmetic in lower precision than the data type? */
	const CCTK_REAL16 (*const input_array_ptr_complex32)[COMPLEX_N_PARTS]
		= ((const CCTK_REAL16 (*)[COMPLEX_N_PARTS]) input_array_ptr)
		  + input_posn;
	LOAD_DATA_COMPLEX32(input_array_ptr_complex32,
			    STRIDE_IJK, part,
			    &data);
	break;
	  }
#endif

default:
CCTK_CVWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): input datatype %d not supported!\n"
"                              (0-origin) input #in=%d"
	   ,
	   (int) input_array_type_codes[in],
	   in);						/*NOTREACHED*/
return UTIL_ERROR_BAD_INPUT;			/*** ERROR RETURN ***/
				/* end of switch (input_array_type_codes[in]) */
				}
			  }

			  }
			/* end of load input array values */
			}


		if (output_arrays[out] != NULL)
		   then {
			/*
			 * interpolate at this point
			 */
			fp result;

			switch	(operation_codes[out])
				{
			#ifdef HAVE_OP_I
			  case 0:
				result = EVALUATE_MOLECULE(&coeffs_I,
							   &data);
				break;
			#endif
			#ifdef HAVE_OP_DX
			  case 1:
				result = inverse_dx
					 * EVALUATE_MOLECULE(&coeffs_dx,
							     &data);
				break;
			#endif
			#ifdef HAVE_OP_DY
			  case 2:
				result = inverse_dy
					 * EVALUATE_MOLECULE(&coeffs_dy,
							     &data);
				break;
			#endif
			#ifdef HAVE_OP_DZ
			  case 3:
				result = inverse_dz
					 * EVALUATE_MOLECULE(&coeffs_dz,
							     &data);
				break;
			#endif
			#ifdef HAVE_OP_DXX
			  case 11:
				result = inverse_dx * inverse_dx
					 * EVALUATE_MOLECULE(&coeffs_dxx,
							     &data);
				break;
			#endif
			#ifdef HAVE_OP_DXY
			  case 12:
			  case 21:
				result = inverse_dx * inverse_dy
					 * EVALUATE_MOLECULE(&coeffs_dxy,
							     &data);
				break;
			#endif
			#ifdef HAVE_OP_DXZ
			  case 13:
			  case 31:
				result = inverse_dx * inverse_dz
					 * EVALUATE_MOLECULE(&coeffs_dxz,
							     &data);
				break;
			#endif
			#ifdef HAVE_OP_DYY
			  case 22:
				result = inverse_dy * inverse_dy
					 * EVALUATE_MOLECULE(&coeffs_dyy,
							     &data);
				break;
			#endif
			#ifdef HAVE_OP_DYZ
			  case 23:
			  case 32:
				result = inverse_dy * inverse_dz
					 * EVALUATE_MOLECULE(&coeffs_dyz,
							     &data);
				break;
			#endif
			#ifdef HAVE_OP_DZZ
			  case 33:
				result = inverse_dz * inverse_dz
					 * EVALUATE_MOLECULE(&coeffs_dzz,
							     &data);
				break;
			#endif
			  default:
CCTK_CVWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): generalized interpolation operation_code %d\n"
"                              is not supported!\n"
"                              (0-origin) output #out=%d"
	   ,
	   (int) operation_codes[out],
	   out);
				return UTIL_ERROR_BAD_INPUT;
							/*** ERROR RETURN ***/
				/* end of switch (operation_codes[out]) */
				}


			/*
			 * store result in output array
			 */
			if (debug >= 10)
			   then {
				if ((pt & (pt-1)) == 0) /* pt is 0 or a power of 2 */
				   then {
					printf("AEILocalInterp:: pt=%d out=%d --> storing result %g\n", pt, out, (double) result);
					fflush(stdout);
					}
				}

			switch	(output_array_type_codes[out])
				{

case CCTK_VARIABLE_REAL:
	  {
	CCTK_REAL *const output_array_ptr_real
		= (CCTK_REAL *) output_arrays[out];
	if (debug >= 10)
	   then {
		if ((pt & (pt-1)) == 0) /* pt is 0 or a power of 2 */
		   then {
			printf("                 result addr is %p\n",
			       (void *) &output_array_ptr_real[pt]);
			printf("                 previous value there was %g\n",
			       output_array_ptr_real[pt]);
			fflush(stdout);
			}
		}
	output_array_ptr_real[pt] = (CCTK_REAL) result;
	break;
	  }

#ifdef HAVE_CCTK_REAL4
case CCTK_VARIABLE_REAL4:
	  {
	CCTK_REAL4 *const output_array_ptr_real4
		= (CCTK_REAL4 *) output_arrays[out];
	output_array_ptr_real4[pt] = (CCTK_REAL4) result;
	break;
	  }
#endif

#ifdef HAVE_CCTK_REAL8
case CCTK_VARIABLE_REAL8:
	  {
	CCTK_REAL8 *const output_array_ptr_real8
		= (CCTK_REAL8 *) output_arrays[out];
	output_array_ptr_real8[pt] = (CCTK_REAL8) result;
	break;
	  }
#endif

#ifdef HAVE_CCTK_REAL16
case CCTK_VARIABLE_REAL16:
	  {
	CCTK_REAL16 *const output_array_ptr_real16
		= (CCTK_REAL16 *) output_arrays[out];
	output_array_ptr_real16[pt] = (CCTK_REAL16) result;
	break;
	  }
#endif

case CCTK_VARIABLE_COMPLEX:
	  {
	CCTK_REAL (*const output_array_ptr_complex)[COMPLEX_N_PARTS]
		= (CCTK_REAL (*)[COMPLEX_N_PARTS]) output_arrays[out];
	output_array_ptr_complex[pt][part] = (CCTK_REAL) result;
	break;
	  }

#ifdef HAVE_CCTK_COMPLEX8
case CCTK_VARIABLE_COMPLEX8:
	  {
	CCTK_REAL4 (*const output_array_ptr_complex8)[COMPLEX_N_PARTS]
		= (CCTK_REAL4 (*)[COMPLEX_N_PARTS]) output_arrays[out];
	output_array_ptr_complex8[pt][part] = (CCTK_REAL4) result;
	break;
	  }
#endif	/* HAVE_CCTK_COMPLEX8 */

#ifdef HAVE_CCTK_COMPLEX16
case CCTK_VARIABLE_COMPLEX16:
	  {
	CCTK_REAL8 (*const output_array_ptr_complex16)[COMPLEX_N_PARTS]
		= (CCTK_REAL8 (*)[COMPLEX_N_PARTS]) output_arrays[out];
	output_array_ptr_complex16[pt][part] = (CCTK_REAL8) result;
	break;
	  }
#endif	/* HAVE_CCTK_COMPLEX16 */

#ifdef HAVE_CCTK_COMPLEX32
case CCTK_VARIABLE_COMPLEX32:
	  {
	CCTK_REAL16 (*const output_array_ptr_complex32)[COMPLEX_N_PARTS]
		= (CCTK_REAL16 (*)[COMPLEX_N_PARTS]) output_arrays[out];
	output_array_ptr_complex32[pt][part] = (CCTK_REAL16) result;
	break;
	  }
#endif	/* HAVE_CCTK_COMPLEX32 */

default:
	CCTK_CVWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): output datatype %d not supported!\n"
"                              (0-origin) output #out=%d"
		   ,
		   (int) output_array_type_codes[out],
		   out);
	return UTIL_ERROR_BAD_INPUT;			/*** ERROR RETURN ***/
				/* end of switch (output type code) */
				}

			/* end of if (output_arrays[out] != NULL) */
			}


		/*
		 * handle querying the Jacobian
		 */
		if ( (Jacobian_info != NULL)
		     && (Jacobian_info->Jacobian_pointer[out] != NULL))
		   then {
			CCTK_REAL *const Jacobian_ptr
			   = Jacobian_info->Jacobian_pointer[out]
			     + Jacobian_info->Jacobian_offset[out]
			     + Jacobian_info->Jacobian_interp_point_stride*pt
			     + Jacobian_info->Jacobian_part_stride*part;

			switch	(operation_codes[out])
				{
			#ifdef HAVE_OP_I
			  case 0:
				STORE_COEFFS(1.0, &coeffs_I,
					     Jacobian_ptr,
					     JACOBIAN_MIJK_STRIDE);
				break;
			#endif
			#ifdef HAVE_OP_DX
			  case 1:
				STORE_COEFFS(inverse_dx, &coeffs_dx,
					     Jacobian_ptr,
					     JACOBIAN_MIJK_STRIDE);
				break;
			#endif
			#ifdef HAVE_OP_DY
			  case 2:
				STORE_COEFFS(inverse_dy, &coeffs_dy,
					     Jacobian_ptr,
					     JACOBIAN_MIJK_STRIDE);
				break;
			#endif
			#ifdef HAVE_OP_DZ
			  case 3:
				STORE_COEFFS(inverse_dz, &coeffs_dz,
					     Jacobian_ptr,
					     JACOBIAN_MIJK_STRIDE);
				break;
			#endif
			#ifdef HAVE_OP_DXX
			  case 11:
				STORE_COEFFS(inverse_dx*inverse_dx, &coeffs_dxx,
					     Jacobian_ptr,
					     JACOBIAN_MIJK_STRIDE);
				break;
			#endif
			#ifdef HAVE_OP_DXY
			  case 12:
			  case 21:
				STORE_COEFFS(inverse_dx*inverse_dy, &coeffs_dxy,
					     Jacobian_ptr,
					     JACOBIAN_MIJK_STRIDE);
				break;
			#endif
			#ifdef HAVE_OP_DXZ
			  case 13:
			  case 31:
				STORE_COEFFS(inverse_dx*inverse_dz, &coeffs_dxz,
					     Jacobian_ptr,
					     JACOBIAN_MIJK_STRIDE);
				break;
			#endif
			#ifdef HAVE_OP_DYY
			  case 22:
				STORE_COEFFS(inverse_dy*inverse_dy, &coeffs_dyy,
					     Jacobian_ptr,
					     JACOBIAN_MIJK_STRIDE);
				break;
			#endif
			#ifdef HAVE_OP_DYZ
			  case 23:
			  case 32:
				STORE_COEFFS(inverse_dy*inverse_dz, &coeffs_dyz,
					     Jacobian_ptr,
					     JACOBIAN_MIJK_STRIDE);
				break;
			#endif
			#ifdef HAVE_OP_DZZ
			  case 33:
				STORE_COEFFS(inverse_dz*inverse_dz, &coeffs_dzz,
					     Jacobian_ptr,
					     JACOBIAN_MIJK_STRIDE);
				break;
			#endif
			  default:
CCTK_CVWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   CCTK_InterpLocalUniform(): generalized interpolation operation_code %d\n"
"                              is not supported!\n"
"                              (0-origin) output #out=%d"
	   ,
	   (int) operation_codes[out],
	   out);
				return UTIL_ERROR_BAD_INPUT;
							/*** ERROR RETURN ***/
				/* end of switch(operation_codes[out])*/
				}
			/* end of Jacobian-query code */
			}

		/* end of  for (part = ...)  loop */
		}
	      }
	  /* end of  for (out = ...)  loop */
	  }
      }
      }
      }
      }

    /* end of  for (pt = ...)  loop */
    }

/*
 * if we're writing logging information, finish this
 */
if (log_fp != NULL)
   then {
	switch	(N_interp_points)
		{
	case 0:
		/* we printed 0 points, but need 4 */
		fprintf(log_fp, "\t0.0\t0.0\t0.0");
		fprintf(log_fp, "\t0.0\t0.0\t0.0");
		/* fall through */
	case 1:
		/* we printed 2 points (0, N_interp_points-1), but need 4 */
		fprintf(log_fp, "\t0.0\t0.0\t0.0");
		/* fall through */
	case 2:
		/* we printed 3 points (0, 1, N_interp_points-1), but need 4 */
		fprintf(log_fp, "\t0.0\t0.0\t0.0");
		/* fall through */
		}

	fprintf(log_fp, "\n");
	fflush(log_fp);
	}

return return_status;
  }
  }
  }
}

/******************************************************************************/

/*@@
  @routine	compute_coeffs_I
  @routine	compute_coeffs_dx
  @routine	compute_coeffs_dy
  @routine	compute_coeffs_dz
  @routine	compute_coeffs_dxx
  @routine	compute_coeffs_dxy
  @routine	compute_coeffs_dxz
  @routine	compute_coeffs_dyy
  @routine	compute_coeffs_dyz
  @routine	compute_coeffs_dzz
  @date		29.Aug.2002
  @author	Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc
	Each of these functions computes the corresponding set of
	interpolation coefficients at the point given by the XYZ
	variables, using the machine-generated experessions in the
	files
		COEFFS_{I,DX,DY,DXX,DXY,DYY,...}_COMPUTE_FILE_NAME

	These functions *must* be static, because this whole file
	is #included (and hence compiled) multiple times.
  @enddesc
  @@*/

#ifdef HAVE_OP_I
  static
    void compute_coeffs_I(FP_XYZ, struct COEFFS_STRUCT *coeffs_I)
  {
  #include COEFFS_I_COMPUTE_FILE_NAME
  }
#endif

#ifdef HAVE_OP_DX
  static
    void compute_coeffs_dx(FP_XYZ, struct COEFFS_STRUCT *coeffs_dx)
  {
  #include COEFFS_DX_COMPUTE_FILE_NAME
  }
#endif

#ifdef HAVE_OP_DY
  static
    void compute_coeffs_dy(FP_XYZ, struct COEFFS_STRUCT *coeffs_dy)
  {
  #include COEFFS_DY_COMPUTE_FILE_NAME
  }
#endif

#ifdef HAVE_OP_DZ
  static
    void compute_coeffs_dz(FP_XYZ, struct COEFFS_STRUCT *coeffs_dz)
  {
  #include COEFFS_DZ_COMPUTE_FILE_NAME
  }
#endif

#ifdef HAVE_OP_DXX
  static
    void compute_coeffs_dxx(FP_XYZ, struct COEFFS_STRUCT *coeffs_dxx)
  {
  #include COEFFS_DXX_COMPUTE_FILE_NAME
  }
#endif

#ifdef HAVE_OP_DXY
  static
    void compute_coeffs_dxy(FP_XYZ, struct COEFFS_STRUCT *coeffs_dxy)
  {
  #include COEFFS_DXY_COMPUTE_FILE_NAME
  }
#endif

#ifdef HAVE_OP_DXZ
  static
    void compute_coeffs_dxz(FP_XYZ, struct COEFFS_STRUCT *coeffs_dxz)
  {
  #include COEFFS_DXZ_COMPUTE_FILE_NAME
  }
#endif

#ifdef HAVE_OP_DYY
  static
    void compute_coeffs_dyy(FP_XYZ, struct COEFFS_STRUCT *coeffs_dyy)
  {
  #include COEFFS_DYY_COMPUTE_FILE_NAME
  }
#endif

#ifdef HAVE_OP_DYZ
  static
    void compute_coeffs_dyz(FP_XYZ, struct COEFFS_STRUCT *coeffs_dyz)
  {
  #include COEFFS_DYZ_COMPUTE_FILE_NAME
  }
#endif

#ifdef HAVE_OP_DZZ
  static
    void compute_coeffs_dzz(FP_XYZ, struct COEFFS_STRUCT *coeffs_dzz)
  {
  #include COEFFS_DZZ_COMPUTE_FILE_NAME
  }
#endif
