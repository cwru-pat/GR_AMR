// patch_interp.cc -- interpolation from a patch
// $Header$

//
// patch_interp::patch_interp
// patch_interp::~patch_interp
// patch_interp::interpolate
//
// patch_interp::verify_Jacobian_sparsity_pattern_ok
// patch_interp::molecule_minmax_ipar_m
// patch_interp::molecule_posn
// patch_interp::Jacobian
// patch_interp::query_interpolator
//

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "../jtutil/util_Table.h"
#include "../jtutil/interpolator/InterpLocalUniform.h"
#include "../cctk.h"

#include "../jtutil/util.hh"
#include "../jtutil/array.hh"
#include "../jtutil/cpm_map.hh"
#include "../jtutil/linear_map.hh"


#include "coords.hh"
#include "grid.hh"
#include "fd_grid.hh"
#include "patch.hh"
#include "patch_edge.hh"
#include "patch_interp.hh"
#include "ghost_zone.hh"

// all the code in this file is inside this namespace
namespace AHFD
	  {

int CCTK_InterpLocalUniform(int N_dims,
                            int operator_handle,
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
  int retval;
  retval =  AEILocalInterp_U_Lagrange_TP(N_dims, param_table_handle,
                          /*** coordinate system ***/
                          coord_origin, coord_delta,
                          /*** interpolation points ***/
                          N_interp_points, interp_coords_type_code,
                          interp_coords,
                          /***** input arrays *****/
                          N_input_arrays, input_array_dims,
                          input_array_type_codes, input_arrays,
                          /***** output arrays *****/
                          N_output_arrays, output_array_type_codes,
                          output_arrays);

  return (retval);
}

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

//
// This function constructs an  patch_interp  object.
//
// Note this requires that the adjacent-side  ghost_zone  objects already
// exist, since they're used in determining the ipar range over which the
// interpolation will be done.
//
// It also requires that the patch's gridfns exist, since we size various
// arrays based on the patch's min/max ghosted gfn.
//
patch_interp::patch_interp(const patch_edge& my_edge_in,
			   int min_iperp_in, int max_iperp_in,
			   const jtutil::array1d<int>& min_parindex_array_in,
			   const jtutil::array1d<int>& max_parindex_array_in,
			   const jtutil::array2d<fp>& interp_par_in,
			   bool ok_to_use_min_par_ghost_zone,
			   bool ok_to_use_max_par_ghost_zone,
			   int interp_handle_in, int interp_par_table_handle_in)
	: my_patch_(my_edge_in.my_patch()),
	  my_edge_(my_edge_in),
	  min_gfn_(my_patch().ghosted_min_gfn()),
	  max_gfn_(my_patch().ghosted_max_gfn()),
	  ok_to_use_min_par_ghost_zone_(ok_to_use_min_par_ghost_zone),
	  ok_to_use_max_par_ghost_zone_(ok_to_use_max_par_ghost_zone),
	  min_iperp_(min_iperp_in), max_iperp_(max_iperp_in),
	  min_ipar_(ok_to_use_min_par_ghost_zone
		    ? my_edge_in.min_ipar_with_corners()
		    : my_edge_in.min_ipar_without_corners()),
	  max_ipar_(ok_to_use_max_par_ghost_zone
		    ? my_edge_in.max_ipar_with_corners()
		    : my_edge_in.max_ipar_without_corners()),
	  min_parindex_array_(min_parindex_array_in),
	  max_parindex_array_(max_parindex_array_in),
	  interp_par_(interp_par_in),
	  interp_handle_(interp_handle_in),
	  interp_par_table_handle_(Util_TableClone(interp_par_table_handle_in)),
	  // note interpolator coordinate origin = min_ipar_
	  // cf. comments in "patch_interp.hh" on gridfn array subscripting
	  gridfn_coord_origin_(my_edge().par_map().fp_of_int(min_ipar_)),
	  gridfn_coord_delta_ (my_edge().par_map().delta_fp()),
	  gridfn_type_codes_      (min_gfn_, max_gfn_),
	  gridfn_data_ptrs_       (min_gfn_, max_gfn_),
	  interp_data_buffer_ptrs_(min_gfn_, max_gfn_) // no comma
{
int status;

// did interpolator-parameter-table cloning work ok?
status = interp_par_table_handle_;
if (status < 0)
   then error_exit(ERROR_EXIT,
"***** patch_interp::patch_interp():\n"
"        can't clone interpolator parmameter table (handle=%d)!\n"
"        error status=%d\n"
,
		   interp_par_table_handle_in,
		   status);					/*NOTREACHED*/

//
// set up the gridfn array stride in the parameter table
// so the interpolator can access the gridfn data arrays
//

const int N_dims = 1;

// stride
const CCTK_INT stride = my_edge().ghosted_par_stride();
status = Util_TableSetIntArray(interp_par_table_handle_,
			       N_dims, &stride,
			       "input_array_strides");
if (status < 0)
   then error_exit(ERROR_EXIT,
"***** patch_interp::patch_interp():\n"
"        can't set gridfn stride in interpolator parmameter table!\n"
"        error status=%d\n"
,
		   status);					/*NOTREACHED*/

//
// set up the gridfn type codes for the interpolator
//
	for (int gfn = min_gfn_ ; gfn <= max_gfn_ ; ++gfn)
	{
	gridfn_type_codes_(gfn) = CCTK_VARIABLE_REAL;	// fp == CCTK_REAL
	}
}

//*****************************************************************************

//
// This function destroys an  patch_interp  object
//
patch_interp::~patch_interp()
{
Util_TableDestroy(interp_par_table_handle_);
}

//*****************************************************************************

//
// This function interpolates the specified range of ghosted gridfns
// at all the coordinates specified when we were constructed.  (It does
// a separate interpolator call for each iperp.)  This function stores
// the interpolated data in  data_buffer(gfn, iperp,parindex) .
//
void patch_interp::interpolate(int ghosted_min_gfn_to_interp,
			       int ghosted_max_gfn_to_interp,
			       jtutil::array3d<fp>& data_buffer)
	const
{
int status;

const int N_dims = 1;
const int N_gridfns = jtutil::how_many_in_range(ghosted_min_gfn_to_interp,
						ghosted_max_gfn_to_interp);
const CCTK_INT N_gridfn_data_points
	= jtutil::how_many_in_range(min_ipar(), max_ipar());
const int interp_coords_type_code = CCTK_VARIABLE_REAL;
const CCTK_INT *gridfn_type_codes_ptr
	= & gridfn_type_codes_(ghosted_min_gfn_to_interp);

//
// do the interpolations at each iperp
//
	for (int iperp = min_iperp() ; iperp <= max_iperp() ; ++iperp)
	{
	//
	// interpolation-point coordinates
	//
	const int min_parindex = min_parindex_array_(iperp);
	const int max_parindex = max_parindex_array_(iperp);
	const CCTK_INT N_interp_points
		= jtutil::how_many_in_range(min_parindex, max_parindex);
	const fp* const interp_coords_ptr = & interp_par_(iperp, min_parindex);
	const void* const interp_coords[N_dims]
		= { static_cast<const void *>(interp_coords_ptr) };

	//
	// pointers to gridfn data to interpolate, and to result buffer
	//
		for (int ghosted_gfn = ghosted_min_gfn_to_interp ;
		     ghosted_gfn <= ghosted_max_gfn_to_interp ;
		     ++ghosted_gfn)
		{
		// set up data pointer to --> (iperp,min_ipar) gridfn
		const int start_irho
			= my_edge().  irho_of_iperp_ipar(iperp, min_ipar());
		const int start_isigma
			= my_edge().isigma_of_iperp_ipar(iperp, min_ipar());
		gridfn_data_ptrs_(ghosted_gfn)
			= static_cast<const void *>(
				& my_patch()
				  .ghosted_gridfn(ghosted_gfn,
						  start_irho, start_isigma)
						   );
		interp_data_buffer_ptrs_(ghosted_gfn)
			= static_cast<void *>(
				& data_buffer(ghosted_gfn, iperp,min_parindex)
					     );
		}
	const void* const* const gridfn_data
		= & gridfn_data_ptrs_(ghosted_min_gfn_to_interp);
	void* const* const interp_buffer
		= & interp_data_buffer_ptrs_(ghosted_min_gfn_to_interp);

	//
	// interpolate for all the parindex values at this iperp
	//
	status = CCTK_InterpLocalUniform(N_dims,
					 interp_handle_,
					 interp_par_table_handle_,
					 &gridfn_coord_origin_,
					 &gridfn_coord_delta_,
					 N_interp_points,
					    interp_coords_type_code,
					    interp_coords,
					 N_gridfns,
					    &N_gridfn_data_points,
					    gridfn_type_codes_ptr,
					    gridfn_data,
					 N_gridfns,
					    gridfn_type_codes_ptr,
					    interp_buffer);

	if (status < 0)
	   then error_exit(ERROR_EXIT,
"***** patch_interp::interpolate():\n"
"        error return %d from interpolator at iperp=%d of [%d,%d]!\n"
"        my_patch()=\"%s\" my_edge()=\"%s\"\n"
,
			   status, iperp, min_iperp(), max_iperp(),
			   my_patch().name(), my_edge().name()); /*NOTREACHED*/
	}
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function verifies that our interpolator has a Jacobian sparsity
// pattern which we grok: at present this means that the molecules are
// fixed-sized hypercubes, with size/shape independent of the interpolation
// coordinates and the floating-point values in the input arrays.
// (N.b. we don't take derivatives in the interpatch interpolation,
// so it doesn't matter of the molecule size/shape depends on any
// operation_codes[] values.)
//
// FIXME: we should also query/check the molecule family!
//
// If the verification is successful, this function is a no-op.  If not,
// this function does an error_exit() and does not return to its caller.
//
void patch_interp::verify_Jacobian_sparsity_pattern_ok()
	const
{
int status1, status2, status3;

//
// query the interpolator
// ... we don't need to set up anything in the parameter table: the
//     Jacobian metainfo is returned by the interpolator on every call,
//     even without any explicit request
//
query_interpolator("verify_Jacobian_sparsity_pattern_ok");

//
// get Jacobian-sparsity query results from parameter table
//
CCTK_INT MSS_is_fn_of_interp_coords, MSS_is_fn_of_input_array_values;
CCTK_INT Jacobian_is_fn_of_input_array_values;
status1 = Util_TableGetInt(interp_par_table_handle_,
			   &MSS_is_fn_of_interp_coords,
			   "MSS_is_fn_of_interp_coords");
status2 = Util_TableGetInt(interp_par_table_handle_,
			   &MSS_is_fn_of_input_array_values,
			   "MSS_is_fn_of_input_array_values");
status3 = Util_TableGetInt(interp_par_table_handle_,
			   &Jacobian_is_fn_of_input_array_values,
			   "Jacobian_is_fn_of_input_array_values");
if ((status1 < 0) || (status2 < 0) || (status3 < 0))
   then error_exit(ERROR_EXIT,
"***** patch_interp::verify_Jacobian_sparsity_pattern_ok():\n"
"        can't get Jacobian sparsity info\n"
"        from interpolator parmameter table!\n"
"        error status1=%d status2=%d status3=%d\n"
,
		   status1, status2, status3);			/*NOTREACHED*/

//
// verify that we grok the Jacobian sparsity pattern
//
if ( MSS_is_fn_of_interp_coords || MSS_is_fn_of_input_array_values
     || Jacobian_is_fn_of_input_array_values )
   then error_exit(ERROR_EXIT,
"***** patch_interp::verify_Jacobian_sparsity_pattern_ok():\n"
"        implementation restriction: we only grok Jacobians with\n"
"        fixed-sized hypercube-shaped molecules, independent of\n"
"        the interpolation coordinates and the floating-point values!\n"
"        MSS_is_fn_of_interp_coords=(int)%d (we only grok 0)\n"
"        MSS_is_fn_of_input_array_values=(int)%d (we only grok 0)\n"
"        Jacobian_is_fn_of_input_array_values=(int)%d (we only grok 0)\n"
,
		   int(MSS_is_fn_of_interp_coords),
		   int(MSS_is_fn_of_input_array_values),
		   int(Jacobian_is_fn_of_input_array_values));	/*NOTREACHED*/
}

//******************************************************************************

//
// This function queries the interpolator to get the [min,max] ipar m
// coordinates of the interpolation molecules.
//
// (This API implicitly assumes that the Jacobian sparsity is one which
// is "ok" as verified by  verify_Jacobian_sparsity_pattern_ok() .)
//
void patch_interp::molecule_minmax_ipar_m(int& min_ipar_m, int& max_ipar_m)
	const
{
const int N_dims = 1;
int status1, status2;

//
// set molecule min/max m query entries in parameter table
//
status1 = Util_TableSetInt(interp_par_table_handle_,
			   0, "molecule_min_m");
status2 = Util_TableSetInt(interp_par_table_handle_,
			   0, "molecule_max_m");
if ((status1 < 0) || (status2 < 0))
   then error_exit(ERROR_EXIT,
"***** patch_interp::molecule_minmax_ipar_m():\n"
"        can't set molecule min/max m queries\n"
"        in interpolator parmameter table!\n"
"        error status1=%d status2=%d\n"
,
		   status1, status2);				/*NOTREACHED*/

//
// query the interpolator
//
query_interpolator("molecule_minmax_ipar_m");

//
// get molecule min/max m query results from parameter table
// ... must go through temp variables for CCTK_INT to int type conversion
//
CCTK_INT molecule_min_m, molecule_max_m;
status1 = Util_TableGetIntArray(interp_par_table_handle_,
				N_dims, &molecule_min_m,
				"molecule_min_m");
status2 = Util_TableGetIntArray(interp_par_table_handle_,
				N_dims, &molecule_max_m,
				"molecule_max_m");
if ((status1 < 0) || (status2 < 0))
   then error_exit(ERROR_EXIT,
"***** patch_interp::molecule_maxmax_ipar_m():\n"
"        can't get molecule min/max m query results\n"
"        from interpolator parmameter table!\n"
"        error status1=%d status2=%d\n"
,
		   status1, status2);				/*NOTREACHED*/
min_ipar_m = molecule_min_m;
max_ipar_m = molecule_max_m;

//
// delete Jacobian-sparsity-pattern query entries from the parameter table
// (so future interpolator calls don't mistakenly redo this query)
//
status1 = Util_TableDeleteKey(interp_par_table_handle_,
			      "molecule_min_m");
status2 = Util_TableDeleteKey(interp_par_table_handle_,
			      "molecule_max_m");
if ((status1 < 0) || (status2 < 0))
   then error_exit(ERROR_EXIT,
"***** patch_interp::verify_Jacobian_sparsity_pattern_ok():\n"
"        can't delete molecule min/max m queries\n"
"        from interpolator parmameter table!\n"
"        error status1=%d status2=%d\n"
,
		   status1, status2);				/*NOTREACHED*/
}

//******************************************************************************

//
// This function queries the interpolator at each iperp to find out the
// molecule ipar positions (which we implicitly assume to be independent
// of ghosted_gfn), and stores these in  posn_buffer(iperp, parindex) .
//
// (This API implicitly assumes that the Jacobian sparsity is one which
// is "ok" as verified by  verify_Jacobian_sparsity_pattern_ok() .)
//
void patch_interp::molecule_posn(jtutil::array2d<CCTK_INT>& posn_buffer)
	const
{
const int N_dims = 1;
int status;

	for (int iperp = min_iperp() ; iperp <= max_iperp() ; ++iperp)
	{
	const int min_parindex = min_parindex_array_(iperp);
	const int max_parindex = max_parindex_array_(iperp);

	// set up the molecule-position query in the parameter table
	CCTK_POINTER molecule_posn_ptrs[N_dims]
	    = { static_cast<CCTK_POINTER>(& posn_buffer(iperp, min_parindex)) };
	status = Util_TableSetPointerArray(interp_par_table_handle_,
					   N_dims, molecule_posn_ptrs,
					   "molecule_positions");
	if (status < 0)
	   then error_exit(ERROR_EXIT,
"***** patch_interp::molecule_posn():\n"
"        can't set molecule position query\n"
"        in interpolator parmameter table at iperp=%d of [%d,%d]!\n"
"        error status=%d\n"
,
			   iperp, min_iperp(), max_iperp(),
			   status);				/*NOTREACHED*/

	// query the interpolator to get the molecule positions
	// for all parindex at this iperp; the interpolator will
	// store the  parindex-min_parindex  values (cf comments
	// on array subscripting at the start of "patch_interp.hh")
	// directly in the  posn_buffer(iperp,*)  array)
	query_interpolator("molecule_positions", iperp);

	// convert the molecule positions from  parindex-min_ipar
	// to  parindex  values (again, cf comments on array subscripting
	// at the start of "patch_interp.hh")
		for (int parindex = min_parindex ;
		     parindex <= max_parindex ;
		     ++parindex)
		{
		posn_buffer(iperp, parindex) += min_ipar();
		}
	}

//
// if we actually did any queries,
// delete molecule-position query entry from the parameter table
// (so future interpolator calls don't mistakenly redo this query)
//
if (jtutil::how_many_in_range(min_iperp(), max_iperp()) > 0)
   then {
	status = Util_TableDeleteKey(interp_par_table_handle_,
				     "molecule_positions");
	if (status < 0)
	   then error_exit(ERROR_EXIT,
"***** patch_interp::molecule_posn():\n"
"        can't delete molecule position query"
"        from interpolator parmameter table!\n"
"        error status=%d\n"
,
			   status);				/*NOTREACHED*/
	}
}

//*****************************************************************************

//
// This function queries the interpolator at each iperp to get the
// Jacobian of the interpolated data with respect to this patch's
// ghosted gridfns,
//	partial interpolate() data_buffer(ghosted_gfn, iperp, parindex)
//	---------------------------------------------------------------
//	    partial ghosted_gridfn(ghosted_gfn, iperp, posn+ipar_m)
//
// This function stores the Jacobian in
//	Jacobian_buffer(iperp, parindex, ipar_m)
// where we implicitly assume the Jacobian to be independent of
// ghosted_gfn[*], and where
//	posn = posn_buffer(iperp, parindex)
// is the molecule position as returned by  molecule_posn() .
//
// [*] We actually do the queries for  min_gfn_  (as computed in our
//     constructor).
//
// This function calls  error_exit() if any interpolation point is out
// of range or any other error is detected.
//
// (This API implicitly assumes that the Jacobian sparsity is one which
// is "ok" as verified by  verify_Jacobian_sparsity_pattern_ok() .)
//
void patch_interp::Jacobian(jtutil::array3d<fp>& Jacobian_buffer)
	const
{
const int N_dims = 1;
const int N_gridfns = 1;

int status, status1, status2;


//
// set Jacobian stride info in parameter table
//
const int Jacobian_interp_point_stride
	= Jacobian_buffer.subscript_stride_j();
status1 = Util_TableSetInt(interp_par_table_handle_,
			   Jacobian_interp_point_stride,
			   "Jacobian_interp_point_stride");
const CCTK_INT Jacobian_m_strides[N_dims]
	= { Jacobian_buffer.subscript_stride_k() };	// = { 1 } in practice
status2 = Util_TableSetIntArray(interp_par_table_handle_,
				N_dims, Jacobian_m_strides,
				"Jacobian_m_strides");
if ((status1 < 0) || (status2 < 0))
   then error_exit(ERROR_EXIT,
"***** patch_interp::Jacobian():\n"
"        can't set Jacobian stride info in interpolator parmameter table!\n"
"        error status1=%d status2=%d\n"
,
		   status1, status2);			/*NOTREACHED*/


//
// query the Jacobians at each iperp
//
	for (int iperp = min_iperp() ; iperp <= max_iperp() ; ++iperp)
	{
	const int min_parindex = min_parindex_array_(iperp);

	//
	// set up the Jacobian query in the parameter table
	//
	CCTK_POINTER const Jacobian_ptrs[N_gridfns]
		= { static_cast<CCTK_POINTER>(
			& Jacobian_buffer(iperp, min_parindex, 0)
					     ) };
	status = Util_TableSetPointerArray(interp_par_table_handle_,
					   N_gridfns, Jacobian_ptrs,
					   "Jacobian_pointer");
	if (status < 0)
	   then error_exit(ERROR_EXIT,
"***** patch_interp::Jacobian():\n"
"        can't set Jacobian pointer in interpolator parmameter table!\n"
"        error status=%d at iperp=%d of [%d,%d]!\n"
,
			   status, iperp, min_iperp(), max_iperp());
								/*NOTREACHED*/

	//
	// query the interpolator to get the Jacobian for this iperp
	//
	query_interpolator("Jacobian", iperp);
	}

//
// delete the query entries from the parameter table
// (so future interpolator calls don't mistakenly redo this query)
//
status = (jtutil::how_many_in_range(min_iperp(), max_iperp()) > 0)
	 ? Util_TableDeleteKey(interp_par_table_handle_,
			       "Jacobian_pointer")
	 : 0;
status1 = Util_TableDeleteKey(interp_par_table_handle_,
			      "Jacobian_interp_point_stride");
status2 = Util_TableDeleteKey(interp_par_table_handle_,
			      "Jacobian_m_strides");
if ((status < 0) || (status1 < 0) || (status2 < 0))
   then error_exit(ERROR_EXIT,
"***** patch_interp::interpolate():\n"
"        can't delete Jacobian query entries\n"
"        from interpolator parmameter table!\n"
"        error status=%d status1=%d status2=%d\n"
,
		   status, status1, status2);			/*NOTREACHED*/
}

//******************************************************************************

//
// This function queries the interpolator with whatever arguments are
// in the parameter table.  It specifies arguments such that no actual
// interpolation is done.
//
// In particular, the following interpolator arguments are set up properly
// based on  iperp :
//	N_dims
//	operator_handle, param_table_handle,
//	coord_origin, coord_delta,
//	N_interp_points, interp_coords_type_code, interp_coords,
//	input_array_dims	# specifies the correct par size of the
//				# patch interpolation region at this iperp
// The following arguments are set to specify a single input array, but
// with a NULL data pointer so no actual data is used:
//	N_input_arrays
//	input_array_type_codes, input_arrays
// The following arguments are set to specify a single output array, but
// with a NULL data pointer so no actual data is stored:
//	N_output_arrays
//	output_array_type_codes, output_arrays
//
// Arguments:
// iperp = (in) Specifies where in the patch interpolation region the
//		interpolator query should be done.
//
// Results:
// If the interpolator returns an error status, this function does an
//  error_exit()  (and does not return to its caller).  Otherwise, this
// function returns the interpolator status code.
//
int patch_interp::query_interpolator(const char function_name[], int iperp)
	const
{
int status;
const int N_dims = 1;

const int min_parindex = min_parindex_array_(iperp);
const int max_parindex = max_parindex_array_(iperp);
const int N_interp_points = jtutil::how_many_in_range(min_parindex,
						      max_parindex);
const int interp_coords_type_code = CCTK_VARIABLE_REAL;
const void* const interp_coords[N_dims]
	= { static_cast<const void *>(& interp_par_(iperp, min_parindex)) };

const int N_input_arrays = 1;
const CCTK_INT input_array_dims[N_dims]
	= { jtutil::how_many_in_range(min_ipar(), max_ipar()) };
const CCTK_INT input_array_type_codes[N_input_arrays] = { CCTK_VARIABLE_REAL };
const void* const input_arrays[N_input_arrays] = { NULL };

const int N_output_arrays = 1;
const CCTK_INT output_array_type_codes[N_output_arrays]
	= { CCTK_VARIABLE_REAL };
void* const output_arrays[N_output_arrays] = { NULL };

status = CCTK_InterpLocalUniform(N_dims,
				 interp_handle_, interp_par_table_handle_,
				 &gridfn_coord_origin_, &gridfn_coord_delta_,
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
if (status < 0)
   then error_exit(ERROR_EXIT,
"***** patch_interp::query_interpolator():\n"
"        on behalf of patch_interp::%s()\n"
"        error return %d from interpolator query at iperp=%d of [%d,%d]!\n"
,
		   function_name,
		   status, iperp, min_iperp(), max_iperp());	/*NOTREACHED*/

return status;
}


            


            
//******************************************************************************
//******************************************************************************
//******************************************************************************

	  }	// namespace AHFinderDirect
