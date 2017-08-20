#ifndef AHFD_PATCH_INTERP_H
#define AHFD_PATCH_INTERP_H
// patch_interp.hh -- interpolation from a patch
// $Header$

//
// prerequisites:
//	<stdio.h>
//	<assert.h>
//	<math.h>
//
//	"cctk.h"
//
//	"stdc.h"
//	"config.hh"
//	"../jtutil/util.hh"
//	"../jtutil/array.hh"
//	"../jtutil/linear_map.hh"
//	"coords.hh"
//	"grid.hh"
//	"fd_grid.hh"
//	"patch.hh"
//

// everything in this file is inside this namespace
namespace AHFD
	  {

            //*****************************************************************************
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
                                        void *const output_arrays[]);

//
// patch_interp - interpolation from a patch
//

//
// A patch_interp object is responsible for interpolating gridfn data
// from its owning patch for use by another patch's ghost_zone object
// (in setting up the gridfn in the other ghost zone).  A patch_interp
// object deals only in its own patch's coordinates; other code elsewhere
// (in practice in interpatch_ghost_zone::) is responsible for translating
// other patch's coordinates into our coordinates.
//

//
// A patch_interp defines a "patch interpolation region", the region of
// its owning patch from which this interpolation will use gridfn data.
//

//
// The way the patch coordnates are constructed, any two adjacent patches
// share a common (perpendicular) coordinate.  Thus we only have to do
// 1-dimensional interpolation here (in the parallel direction).  In
// other words, for each iperp we interpolate in par.
//
// In general we interpolate each gridfn at a number of distinct par
// for each iperp; the integer "parindex" indexes these points.  We
// attach no particular semantics to parindex, and it need not be
// 0-origin or have the same range for each iperp.  [In practice,
// parindex will be the other patch's ipar coordinate.]  However,
// we assume that the range of parindex is roughly similar for each
// iperp, so it's ok to use (iperp,parindex) as a 2-D rectangular
// index space.
//
// For example, we might interpolate at the points
//            ipar ipar ipar ipar ipar ipar ipar ipar ipar
//              1    2    3    4    5    6    7    8    9
// iperp=10           (2a)   (3b)   (4c)
// iperp=11          (2d)   (3e)  (4f)   (5g)
// where the (2a)-(5g) are the interpolation points, with 2-5 being the
// parindex values and a-g being unique identifiers used in our description
// below.  In terms of our member data, this interpolation region would
// be described by
//	[min,max]_iperp_=[10,11]
//	[min,max]_ipar_=[1,9]
//	[min,max]_parindex_array_(10)=[2,5]
//	[min,max]_parindex_array_(11)=[2,6]
//	interp_par_(10,2) = x[a]
//	interp_par_(10,3) = x[b]
//	interp_par_(10,4) = x[c]
//	interp_par_(11,2) = x[d]
//	interp_par_(11,3) = x[e]
//	interp_par_(11,4) = x[f]
//	interp_par_(11,5) = x[g]
//

//
// We use the Cactus local interpolator CCTK_InterpLocalUniform()
// to do the interpolation.  To minimize interpolator overheads, we
// interpolate all the gridfns at each iperp in a single interpolator
// call.  [Different iperp values involve different sets of (1-D)
// gridfn data, and so inherently require distinct interpolator calls.]
//
// Setting up the array subscripting for the interpolator to access
// the gridfn data is a bit tricky:  The interpolator accesses the
// gridfn data using the generic (1-D) subscripting expression
//	data[offset + i*stride]
// where  i  is the data array index.  However, we'd rather not use
//  offset , because it has to be supplied in the parameter table as
// an array subscripted by  gfn , and so would require changing the
// parameter table for each call on  interpolate()  (with potentially
// different numbers of gridfns being interpolated).  Instead, at each
//  iperp  we use  i = ipar-min_ipar , so the default  offset=0  makes
// the subscripting expression zero for  ipar = min_ipar .  This also
// makes the interpolator's  min_i = 0  and  max_i  be  dims-1  (both
// the defaults), so those also don't have to be set in the parameter
// table either.  We set the interpolator's data coordinate origin to
// the  par  coordinate for  min_ipar , so it correctly maps  i --> par .
// With this strategy we can share the interpolator parameter table
// across all the  iperp  values, and we don't need to modify the
// parameter table at all after the initial setup in our constructor.
// However, we do have to adjust the molecule positions in
//  patch_interp::molecule_posn() , since the interpolator will return
//  i  values, while  molecule_posn()  needs  ipar  values.
//

class	patch_interp
	{
public:
	// to which patch/edge do we belong?
	const patch& my_patch() const { return my_patch_; }
	const patch_edge& my_edge() const { return my_edge_; }


public:
	//
	// ***** main client interface *****
	//
	// interpolate specified range of ghosted gridfns
	// at all the coordinates specified when we were constructed,
	// store interpolated data in
	//	data_buffer(ghosted_gfn, iperp, parindex)
	void interpolate(int ghosted_min_gfn_to_interp,
			 int ghosted_max_gfn_to_interp,
			 jtutil::array3d<fp>& data_buffer)
		const;



public:
	//
	// ***** Jacobian of interpolate() *****
	//

	// verify (no-op if ok, error_exit() if not) that interpolator
	// has a Jacobian sparsity pattern which we grok: at present this
	// means molecules are fixed-sized hypercubes, with size/shape
	// independent of interpolation coordinates and the floating-point
	// values in the input arrays
	void verify_Jacobian_sparsity_pattern_ok() const;

	//
	// The API for the remaining Jacobian functions implicitly
	// assumes that the Jacobian sparsity pattern is "ok" as
	// verified by  verify_Jacobian_sparsity_pattern_ok() ,
	// and in particular that  [min,max]_ipar_m  are independent
	// of iperp and parindex.
	//

	// get [min,max] ipar m coordinates of interpolation molecules
	void molecule_minmax_ipar_m(int& min_ipar_m, int& max_ipar_m) const;

	// get interpolation molecule ipar positions in
	//  molecule_posn_buffer(iperp, parindex)
	// ... array type is CCTK_INT so we can pass by reference
	//     to interpolator
	void molecule_posn(jtutil::array2d<CCTK_INT>& posn_buffer) const;

	// get Jacobian of interpolated data with respect to this patch's
	// ghosted gridfns,
	//	partial interpolate() data_buffer(ghosted_gfn, iperp, parindex)
	//	---------------------------------------------------------------
	//	    partial ghosted_gridfn(ghosted_gfn, iperp, posn+ipar_m)
	// store Jacobian in
	//	Jacobian_buffer(iperp, parindex, ipar_m)
	// where we implicitly assume the Jacobian to be independent of
	// ghosted_gfn, and where
	//	posn = posn_buffer(iperp, parindex)
	// as returned by  molecule_posn()
	void Jacobian(jtutil::array3d<fp>& Jacobian_buffer) const;

	//
	// ***** internal functions *****
	//
private:
	// [min,max] iperp for interpolation and gridfn data
	int min_iperp() const { return min_iperp_; }
	int max_iperp() const { return max_iperp_; }

	// min/max (iperp,ipar) of the gridfn data to use for interpolation
	int min_ipar() const { return min_ipar_; }
	int max_ipar() const { return max_ipar_; }

	// query the interpolator via the parameter table
	// on behalf of the specified function (name used for error msgs only)
	// ... error_exit() if error return from interpolator,
	//     otherwise return interpolator status code
	int query_interpolator(const char function_name[], int iperp)
		const;
	// ... use default iperp for queries that don't care
	int query_interpolator(const char function_name[])
		const
		{ return query_interpolator(function_name, min_iperp_); }


	//
	// ***** constructor, destructor, et al *****
	//
public:
	//
	// Constructor arguments:
	// my_edge_in = Identifies the patch/edge to which this
	//		interpolation region is to belong.
	// [min,max]_iperp_in = The range of iperp for this interpolation
	//			region
	// [min,max]_parindex_array_in(iperp)
	//	= [min,max] range of parindex actually used at each iperp.
	//	  We keep references to these arrays, so they should have
	//	  lifetimes at last as long as that of this object.
	// interp_par_in(iperp,parindex)
	//	= Gives the par coordinates at which we will interpolate;
	//	  array entries outside the range [min,max]_parindex_in
	//	  are unused.  We keep a reference to this array, so it
	//	  should have a lifetime at last as long as that of this
	//	  object.
	// ok_to_use_[min,max]_par_ghost_zone
	//	= Boolean flags saying whether or not we should use gridfn
	//	  data from the [min,max]_par ghost zones in the interpolation.
	// interp_handle_in = Cactus handle to the interpatch interpolation
	//		      operator.
	// interp_par_table_handle_in
	//	= Cactus handle to a Cactus key/value table giving
	//	  parameters (eg order) for the interpatch interpolation
	//	  operator.  This class internally clones this table and
	//	  modifies the clone, so the original table is not modified
	//	  by any actions of this class.
	//
	// This constructor requires that this patch's gridfns already
	// exist, since we size various arrays based on the patch's min/max
	// ghosted gfn.
	//
	patch_interp(const patch_edge& my_edge_in,
		     int min_iperp_in, int max_iperp_in,
		     const jtutil::array1d<int>& min_parindex_array_in,
		     const jtutil::array1d<int>& max_parindex_array_in,
		     const jtutil::array2d<fp>& interp_par_in,
		     bool ok_to_use_min_par_ghost_zone,
		     bool ok_to_use_max_par_ghost_zone,
		     int interp_handle_in, int interp_par_table_handle_in);

	~patch_interp();

private:
        // we forbid copying and passing by value
        // by declaring the copy constructor and assignment operator
        // private, but never defining them
	patch_interp(const patch_interp& rhs);
	patch_interp& operator=(const patch_interp& rhs);


	//
	// ***** data members *****
	//
private:
	const patch& my_patch_;
	const patch_edge& my_edge_;

	// range of gfn we can handle
	// (any given interpolate() call may specify a subrange)
	const int min_gfn_, max_gfn_;

	// these are strictly speaking redundant
	// but we keep them for use in debugging
	bool ok_to_use_min_par_ghost_zone_, ok_to_use_max_par_ghost_zone_;

	// patch interpolation region,
	// i.e. range of (iperp,ipar) in this patch from which
	// we will use gridfn data in interpolation
	const int min_iperp_, max_iperp_;
	const int min_ipar_, max_ipar_;

	// [min,max] parindex at each iperp
	// ... these are references to arrays passed in to our constructor
	//     ==> we do *not* own them!
	// ... indices are (iperp)
	const jtutil::array1d<int>& min_parindex_array_;
	const jtutil::array1d<int>& max_parindex_array_;

	// interp_par_(iperp,parindex)
	//	= Gives the par coordinates at which we will interpolate;
	//	  array entries outside the range [min,max]_parindex_in
	//	  are unused (n.b. this interface implicitly takes the
	//	  par coordinates to be independent of ghosted_gfn).
	// ... this is a reference to an array passed in to our constructor
	//     ==> we do *not* own this!
	const jtutil::array2d<fp>& interp_par_;	// indices (iperp,parindex)

	// Cactus handle to the interpolation operator
	int interp_handle_;

	// Cactus handle to our private Cactus key/value table
	// giving parameters for the interpolation operator
	// ... this starts out as a copy of the passed-in table,
	//     then gets extra stuff added to it specific to this
	//     interpolation region; it's shared across all iperp
	// ... we own this table
	const int interp_par_table_handle_;

	// (par) origin and delta values of the gridfn data
	const fp gridfn_coord_origin_, gridfn_coord_delta_;

	// array of gridfn type codes for interpolator
	// ... must be CCTK_INT so we can pass by reference to interpolator
	// ... values set in ctor body, const thereafter
	// ... index is (gfn)
	mutable jtutil::array1d<CCTK_INT> gridfn_type_codes_;

	// --> start of gridfn data to use for interpolation
	//     (reset for each iperp)
	// ... we do *not* own the pointed-to data!
	// ... index is (gfn)
	mutable jtutil::array1d<const void*> gridfn_data_ptrs_;

	// --> start of interpolation data buffer for each gridfn
	//     (reset for each iperp)
	// ... we do *not* own the pointed-to data!
	// ... index is (gfn)
	mutable jtutil::array1d<void*> interp_data_buffer_ptrs_;
	};

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif
