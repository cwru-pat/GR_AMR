#ifndef AHFD_PATCH_PATCH_SYSTEM_H
#define AHFD_PATCH_PATCH_SYSTEM_H
// patch_system.hh -- describes the (an) entire system of interlinked patches
// $Header$

//
// patch_system - describes a system of interlinked patches
//

//
// prerequisites:
//	<stdio.h>
//	<assert.h>
//	<math.h>
//	"stdc.h"
//	"config.hh"
//	"../jtutil/util.hh"
//	"../jtutil/array.hh"
//	"../jtutil/linear_map.hh"
//	"coords.hh"
//	"grid.hh"
//	"fd_grid.hh"
//	"patch.hh"
//	"patch_edge.hh"
//	"ghost_zone.hh"
//	"patch_interp.hh"
//	"patch_info.hh"
//

// everything in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

//
// A  patch_system  object describes a system of interlinked patches.
//
// Its  const  qualifiers refer (only) to the gridfn data.  Notably, this
// means that  synchronize()  is a non-const function (it modifies gridfn
// data), while  synchronize_Jacobian()  et al are const functions (they
// don't modify gridfn data) even though they may update other internal
// state in the  patch_system  object and its subobjects.
//

class	patch_system
	{
	//
	// ***** static data & functions describing patch systems *****
	//
public:
	// what patch-system type are supported?
	// (see "patch_system_info.hh" for detailed descriptions of these)
	enum patch_system_type {
			       patch_system__full_sphere,
			       patch_system__plus_z_hemisphere,
			       patch_system__plus_xy_quadrant_mirrored,
			       patch_system__plus_xy_quadrant_rotating,
			       patch_system__plus_xz_quadrant_mirrored,
			       patch_system__plus_xz_quadrant_rotating,
			       patch_system__plus_xyz_octant_mirrored,
			       patch_system__plus_xyz_octant_rotating
			       };

	// maximum number of patches in any patch-system type
	static
	  const int max_N_patches = 6;

	// decode patch system type into N_patches
	static
	  int N_patches_of_type(enum patch_system_type type_in);

	// patch system type <--> human-readable character-string name
	static
	  const char* name_of_type(enum patch_system_type type_in);
	static
	  enum patch_system_type type_of_name(const char* name_in);


	//
	// ***** coordinates *****
	//
public:
	#ifdef NOT_USED
	// global (x,y,z) --> local (x,y,z)
	fp local_x_of_global_x(fp global_x) const
		{ return global_coords_.local_x_of_global_x(global_x); }
	fp local_y_of_global_y(fp global_y) const
		{ return global_coords_.local_y_of_global_y(global_y); }
	fp local_z_of_global_z(fp global_z) const
		{ return global_coords_.local_z_of_global_z(global_z); }
	#endif	/* NOT_USED */

	#ifdef NOT_USED
	// local (x,y,z) --> global (x,y,z)
	fp global_x_of_local_x(fp local_x) const
		{ return global_coords_.global_x_of_local_x(local_x); }
	fp global_y_of_local_y(fp local_y) const
		{ return global_coords_.global_y_of_local_y(local_y); }
	fp global_z_of_local_z(fp local_z) const
		{ return global_coords_.global_z_of_local_z(local_z); }
	#endif	/* NOT_USED */

	// get global (x,y,z) coordinates of local origin point
	fp origin_x() const { return global_coords_.origin_x(); }
	fp origin_y() const { return global_coords_.origin_y(); }
	fp origin_z() const { return global_coords_.origin_z(); }

	// set global (x,y,z) coordinates of local origin point
        void origin_x(const fp x) { global_coords_.origin_x(x); }
	void origin_y(const fp y) { global_coords_.origin_y(y); }
	void origin_z(const fp z) { global_coords_.origin_z(z); }

 
	//
	// ***** meta-info about the entire patch system *****
	//
public:
	// patch-system type
	enum patch_system_type type() const { return type_; }

	// total number of patches
	int N_patches() const { return N_patches_; }

	// get patches by patch number
	const patch& ith_patch(int pn) const
		{ return * all_patches_[pn]; }
	patch& ith_patch(int pn)
		{ return * all_patches_[pn]; }

	// find a patch by +/- xyz "ctype"
	// FIXME: the present implementation of this function is quite slow
	const patch& plus_or_minus_xyz_patch(bool is_plus, char ctype)
		const;

	// find a patch by name, return patch number; error_exit() if not found
	int patch_number_of_name(const char* name) const;

	// total number of grid points
	int N_grid_points() const { return N_grid_points_; }
	int ghosted_N_grid_points() const { return ghosted_N_grid_points_; }
	int max_N_additional_points() const { return max_N_additional_points_; }
	int N_additional_points() const { return N_additional_points_; }
        int N_additional_points(const int n) { return N_additional_points_=n; }


	//
	// ***** meta-info about gridfns *****
	//
public:
	int min_gfn() const { return ith_patch(0).min_gfn(); }
	int max_gfn() const { return ith_patch(0).max_gfn(); }
	int N_gridfns() const { return ith_patch(0).N_gridfns(); }
	bool is_valid_gfn(int gfn) const
		{ return ith_patch(0).is_valid_gfn(gfn); }
	int ghosted_min_gfn() const { return ith_patch(0).ghosted_min_gfn(); }
	int ghosted_max_gfn() const { return ith_patch(0).ghosted_max_gfn(); }
	int ghosted_N_gridfns() const
		{ return ith_patch(0).ghosted_N_gridfns(); }
	bool is_valid_ghosted_gfn(int ghosted_gfn) const
		{ return ith_patch(0).is_valid_ghosted_gfn(ghosted_gfn); }


	//
	// ***** synchronize() and its Jacobian *****
	//
public:
	// "synchronize" all ghost zones of all patches,
	// i.e. update the ghost-zone values of the specified gridfns
	// via the appropriate sequence of symmetry operations
	// and interpatch interpolations
	void synchronize(int ghosted_min_gfn_to_sync,
			 int ghosted_max_gfn_to_sync);

	// ... do this for all ghosted gridfns
	void synchronize()
		{
		synchronize(ghosted_min_gfn(),
			    ghosted_max_gfn());
		}


	//
	// do any precomputation necessary to compute Jacobian of
	//  synchronize() , taking into account synchronize()'s
	// full 3-phase algorithm
	//
	void compute_synchronize_Jacobian(int ghosted_min_gfn_to_sync,
					  int ghosted_max_gfn_to_sync)
		const;

	// ... do this for all ghosted gridfns
	void compute_synchronize_Jacobian()
		const
		{
		compute_synchronize_Jacobian(ghosted_min_gfn(),
					     ghosted_max_gfn());
		}

	//
	// The following functions access the Jacobian computed by
	//  compute_synchronize_Jacobian() .  Note this API is rather
	// different than that of ghost_zone::comute_Jacobian()  et al:
	// here we must take into account synchronize()'s full 3-phase
	// algorithm, and this may lead to a more general Jacobian
	// structure.
	//
	// This API still implicitly assumes that the Jacobian is
	// independent of  ghosted_gfn , and that the set of y points
	// (with nonzero Jacobian values) in a single row of the Jacobian
	// matrix (i.e. the set of points on which a single ghost-zone
	// point depends),
	// - lies entirely within a single y patch
	// - has a single yiperp value
	// - have a contiguous interval of yipar; we parameterize this
	//   interval as  yipar = posn+m
	//

	// what are the global min/max  m  over all ghost zone points?
	// (this is useful for sizing the buffer for synchronize_Jacobian())
	void synchronize_Jacobian_global_minmax_ym(int& min_ym, int& max_ym)
		const;

	// compute a single row of the Jacobian:
	// - return value is edge to which y point belongs
	//   (caller can get patch from this edge)
	// - store y_iperp and y_posn and min/max ym in named arguments
	// - stores the Jacobian elements
	//	partial synchronize() gridfn(ghosted_gfn, px, x_iperp, x_ipar)
	//	-------------------------------------------------------------
	//	     partial gridfn(ghosted_gfn, py, y_iperp, y_posn+ym)
	//   (taking into account synchronize()'s full 3-phase algorithm)
	//   in the caller-supplied buffer
	//	Jacobian_buffer(ym)
	//   for each  ym  in the min/max ym range
	const patch_edge&
	  synchronize_Jacobian(const ghost_zone& xgz,
			       int x_iperp, int x_ipar,
			       int& y_iperp,
			       int& y_posn, int& min_ym, int& max_ym,
			       jtutil::array1d<fp>& Jacobian_buffer)
		const;

	// helper functions for synchronize_Jacobian():
private:
	// "fold" (part of) a Jacobian row
	// to take a symmetry operation into acount
	// e_Jac = edge which the Jacobian lies along
	// e_fold = edge about which to fold
	// [min,max]_m = range of m in the Jacobian
	// [min,max]_fold_m = range of m to fold
	//		      (must be a subrange of {min,max}_m)
	void fold_Jacobian(const patch_edge& e_Jac, const patch_edge& e_fold,
			   int iperp,
			   int posn, int min_m, int max_m,
			   int min_fold_m, int max_fold_m,
			   jtutil::array1d<fp>& Jacobian_buffer)
		const;

	// compute the Jacobian of ghost zone's synchronize()
	// *without* taking into account 3-phase algorithm
	const patch_edge&
	  ghost_zone_Jacobian(const ghost_zone& xgz,
			      int x_iperp, int x_ipar,
			      int& y_iperp,
			      int& y_posn, int& min_ym, int& max_ym,
			      jtutil::array1d<fp>& Jacobian_buffer)
		const;


	//
	// ***** gridfn operations *****
	//
public:
	// dst = a
	void         set_gridfn_to_constant(fp a, int         dst_gfn);
	void set_ghosted_gridfn_to_constant(fp a, int ghosted_dst_gfn);

	// dst = src
	void         gridfn_copy(int         src_gfn, int         dst_gfn);
	void ghosted_gridfn_copy(int ghosted_src_gfn, int ghosted_dst_gfn);
	void ghosted_gridfn_copy_to_nominal(int ghosted_src_gfn, int dst_gfn);

	// dst += delta
	void         add_to_gridfn(fp delta, int         dst_gfn);
	void add_to_ghosted_gridfn(fp delta, int ghosted_dst_gfn);

	// dst *= factor
	void scale_ghosted_gridfn(fp factor, int ghosted_dst_gfn);

	// compute norms of gridfn (only over nominal grid)
	void         gridfn_norms(int         src_gfn, jtutil::norm<fp>& norms)
		const;
	void ghosted_gridfn_norms(int ghosted_src_gfn, jtutil::norm<fp>& norms)
		const;


	//
	// ***** testing (x,y,z) point position versus a surface *****
	//

	// find patch containing (ray from origin to) given local (x,y,z)
	// ... if there are multiple patches containing the position,
	//     we return the one which would still contain it if patches
	//     didn't overlap; if multiple patches satisfy this criterion
	//     then it's arbitrary which one we return
	// ... if no patch contains the position (for a non--full-sphere
	//     patch system), or the position is at the origin, then
	//     we return a NULL pointer
	const patch* patch_containing_local_xyz(fp x, fp y, fp z)
		const;

	// radius of surface in direction of an (x,y,z) point,
	// taking into account any patch-system symmetries;
	// or dummy value 1.0 if point is identical to local origin
	fp radius_in_local_xyz_direction(int ghosted_radius_gfn,
					 fp x, fp y, fp z)
		const;

	void radii_in_local_xyz_directions(int ghosted_radius_gfn,
                                           int npoints,
                                           fp const * xp,
                                           fp const * yp,
                                           fp const * zp,
                                           fp * radii)
		const;


	//
	// ***** line/surface operations *****
	//

	// compute the circumference of a surface in the {xy, xz, yz} plane
	// ... note this is the full circumference all around the sphere,
	//     even if the patch system only covers a proper subset of this
	// ... the implementation assumes adjacent patches are butt-joined
	// ... plane must be one of "xy", "xz", or "yz"
	fp circumference(const char plane[],
			 int ghosted_radius_gfn,
			 int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
				       int g_yy_gfn, int g_yz_gfn,
						     int g_zz_gfn,
			 enum patch::integration_method method)
		const;

	// compute the surface integral of a gridfn over the 2-sphere
	//	$\int f(\rho,\sigma) \, dA$
	//		= \int f(\rho,\sigma) \sqrt{|J|} \, d\rho \, d\sigma$
	// where $J$ is the Jacobian of $(x,y,z)$ with respect to $(rho,sigma)
	// ... integration method selected by  method  argument
	// ... src gridfn may be either nominal-grid or ghosted-grid
	// ... Boolean flags  src_gfn_is_even_across_{xy,xz,yz}_planes
	//     specify whether the gridfn to be integrated is even (true)
	//     or odd (false) across the corresponding planes.  Only the
	//     flags corresponding to boundaries of the patch system are
	//     used.  For example, for a  plus_z_hemisphere  patch system,
	//     only the  src_gfn_is_even_across_xy_plane  flag is used.
	// ... note integral is over the full 2-sphere,
	//     even if the patch system only covers a proper subset of this
	// ... the implementation assumes adjacent patches are butt-joined
	fp integrate_gridfn(int unknown_src_gfn,
			    bool src_gfn_is_even_across_xy_plane,
			    bool src_gfn_is_even_across_xz_plane,
			    bool src_gfn_is_even_across_yz_plane,
			    int ghosted_radius_gfn,
			    int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
					  int g_yy_gfn, int g_yz_gfn,
							int g_zz_gfn,
			    enum patch::integration_method method)
		const;

	fp integrate_gridpoint(int unknown_src_gfn,
                               int ghosted_radius_gfn,
                               int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
                                             int g_yy_gfn, int g_yz_gfn,
                                                           int g_zz_gfn,
                               enum patch::integration_method method,
                               int pn, int irho, int isigma)
		const;

	fp integrate_correction(bool src_gfn_is_even_across_xy_plane,
                                bool src_gfn_is_even_across_xz_plane,
                                bool src_gfn_is_even_across_yz_plane)
		const;



	//
	// ***** I/O *****
	//
public:
	// print to a named file (newly (re)created)
	// output format is
	//	dpx	dpy	gridfn
	void print_gridfn(int gfn, const char output_file_name[]) const
		{
		print_unknown_gridfn(false, gfn,
				     false, false, 0,
				     output_file_name, false);
		}
	void print_ghosted_gridfn(int ghosted_gfn,
				  const char output_file_name[],
				  bool want_ghost_zones = true)
		const
		{
		print_unknown_gridfn(true, ghosted_gfn,
				     false, false, 0,
				     output_file_name, want_ghost_zones);
		}

	// print to a named file (newly (re)created)
	// output format is
	//	dpx	dpy	gridfn   global_x   global_y   global_z
	// where global_[xyz} are derived from the angular position
	// and a specified (unknown-grid) radius gridfn
	void print_gridfn_with_xyz
		(int gfn,
		 bool radius_is_ghosted_flag, int unknown_radius_gfn,
		 const char output_file_name[])
		const
		{
		print_unknown_gridfn(false, gfn,
				     true, radius_is_ghosted_flag,
					   unknown_radius_gfn,
				     output_file_name, false);
		}
	void print_ghosted_gridfn_with_xyz
		(int ghosted_gfn,
		 bool radius_is_ghosted_flag, int unknown_radius_gfn,
		 const char output_file_name[],
		 bool want_ghost_zones = true)
		const
		{
		print_unknown_gridfn(true, ghosted_gfn,
				     true, radius_is_ghosted_flag,
					   unknown_radius_gfn,
				     output_file_name, want_ghost_zones);
		}

public:
	// read from a named file
	void read_gridfn(int gfn, const char input_file_name[])
		{ read_unknown_gridfn(false, gfn, input_file_name, false); }
	void read_ghosted_gridfn(int ghosted_gfn,
				 const char input_file_name[],
				 bool want_ghost_zones = true)
		{
		read_unknown_gridfn(true, ghosted_gfn,
				    input_file_name, want_ghost_zones);
		}

private:
	// ... internal worker functions
	void print_unknown_gridfn
		(bool ghosted_flag, int unknown_gfn,
		 bool print_xyz_flag, bool radius_is_ghosted_flag,
				      int unknown_radius_gfn,
		 const char output_file_name[], bool want_ghost_zones)
		const;
	void read_unknown_gridfn(bool ghosted_flag, int unknown_gfn,
				 const char input_file_name[],
				 bool want_ghost_zones);

public:
	// output to a named file
	// output format is HDF5
        void output_gridfn(int gfn,
                           const char gfn_name[],
                           const char output_file_name[]) const
		{
		output_unknown_gridfn(false, gfn,
                                      gfn_name,
                                      false, false, 0,
                                      output_file_name, false);
		}
	void output_ghosted_gridfn(int ghosted_gfn,
                                   const char gfn_name[], 
                                   const char output_file_name[],
                                   bool want_ghost_zones = true)
		const
		{
		output_unknown_gridfn(true, ghosted_gfn,
                                      gfn_name,
                                      false, false, 0,
                                      output_file_name, want_ghost_zones);
		}

	// output to a named file (newly (re)created)
	// output format is
	//	dpx	dpy	gridfn   global_x   global_y   global_z
	// where global_[xyz} are derived from the angular position
	// and a specified (unknown-grid) radius gridfn
	void output_gridfn_with_xyz
		(int gfn,
                 const char gfn_name[],
		 bool radius_is_ghosted_flag, int unknown_radius_gfn,
		 const char output_file_name[])
		const
		{
		output_unknown_gridfn(false, gfn,
                                      gfn_name,
                                      true, radius_is_ghosted_flag,
                                            unknown_radius_gfn,
                                      output_file_name, false);
		}
	void output_ghosted_gridfn_with_xyz
		(int ghosted_gfn,
                 const char gfn_name[], 
		 bool radius_is_ghosted_flag, int unknown_radius_gfn,
		 const char output_file_name[],
		 bool want_ghost_zones = true)
		const
		{
		output_unknown_gridfn(true, ghosted_gfn,
                                      gfn_name,
                                      true, radius_is_ghosted_flag,
                                            unknown_radius_gfn,
                                      output_file_name, want_ghost_zones);
		}

private:
	// ... internal worker functions
	void output_unknown_gridfn
		(bool ghosted_flag, int unknown_gfn,
                 const char gfn_name[], 
		 bool print_xyz_flag, bool radius_is_ghosted_flag,
				      int unknown_radius_gfn,
		 const char output_file_name[], bool want_ghost_zones)
		const;



	//
	// ***** access to gridfns as 1-D arrays *****
	//
	// ... n.b. this interface implicitly assumes that gridfn data
	//     arrays are contiguous across patches; this is ensured by
	//     setup_gridfn_storage() (called by our constructor)
	//
public:
	// convert (patch,irho,isigma) <--> 1-D 0-origin grid point number (gpn)
	int gpn_of_patch_irho_isigma(const patch& p, int irho, int isigma)
		const
		{
		return starting_gpn_[p.patch_number()]
		       + p.gpn_of_irho_isigma(irho,isigma);
		}
	int ghosted_gpn_of_patch_irho_isigma(const patch& p,
					     int irho, int isigma)
		const
		{
		return ghosted_starting_gpn_[p.patch_number()]
		       + p.ghosted_gpn_of_irho_isigma(irho,isigma);
		}
	// ... n.b. we return patch as a reference via the function result;
	//     an alternative would be to have a patch*& argument
	const patch&
	  patch_irho_isigma_of_gpn(int gpn, int& irho, int& isigma)
		const;
	const patch&
	  ghosted_patch_irho_isigma_of_gpn(int gpn, int& irho, int& isigma)
		const;

	// access actual gridfn data arrays
	// (low-level, dangerous, use with caution)
	const fp* gridfn_data(int gfn) const
		{ return ith_patch(0).gridfn_data_array(gfn); }
	      fp* gridfn_data(int gfn)
		{ return ith_patch(0).gridfn_data_array(gfn); }
	const fp* ghosted_gridfn_data(int ghosted_gfn) const
		{ return ith_patch(0).ghosted_gridfn_data_array(ghosted_gfn); }
	      fp* ghosted_gridfn_data(int ghosted_gfn)
		{ return ith_patch(0).ghosted_gridfn_data_array(ghosted_gfn); }


	//
	// ***** constructor, destructor *****
	//
	// This constructor doesn't support the full generality of the
	// patch data structures (which would, eg, allow ghost_zone_width
	// and patch_extend_width and the interpolator parameters to vary
	// from ghost zone to ghost zone, and the grid spacings to vary
	// from patch to patch.  But in practice we'd probably never
	// use that generality...
	//
public:
	patch_system(fp origin_x_in, fp origin_y_in, fp origin_z_in,
		     enum patch_system_type type_in,
		     int ghost_zone_width, int patch_overlap_width,
		     int N_zones_per_right_angle,
		     int min_gfn_in, int max_gfn_in,
		     int ghosted_min_gfn_in, int ghosted_max_gfn_in,
		     int ip_interp_handle_in, int ip_interp_par_table_handle_in,
		     int surface_interp_handle_in,
		     int surface_interp_par_table_handle_in,
		     bool print_summary_msg_flag, bool print_detailed_msg_flag);
	~patch_system();


	//
	// ***** helper functions for constructor *****
	//
private:
	// construct patches as described by patch_info[] array,
	// and link them into the patch system
	// does *NOT* create ghost zones
	// does *NOT* set up gridfns
	void create_patches(const struct patch_info patch_info_in[],
			    int ghost_zone_width, int patch_extend_width,
			    int N_zones_per_right_angle,
			    bool print_msg_flag);

	// setup all gridfns with contiguous-across-patches storage
	void setup_gridfn_storage
		(int min_gfn_in, int max_gfn_in,
		 int ghosted_min_gfn_in, int ghosted_max_gfn_in,
		 bool print_msg_flag);

	// setup (create/interlink) all ghost zones
	void setup_ghost_zones__full_sphere
		(int patch_overlap_width,
		 int ip_interp_handle, int ip_interp_par_table_handle,
		 bool print_msg_flag);
	void setup_ghost_zones__plus_z_hemisphere
		(int patch_overlap_width,
		 int ip_interp_handle, int ip_interp_par_table_handle,
		 bool print_msg_flag);
	void setup_ghost_zones__plus_xy_quadrant_mirrored
		(int patch_overlap_width,
		 int ip_interp_handle, int ip_interp_par_table_handle,
		 bool print_msg_flag);
	void setup_ghost_zones__plus_xy_quadrant_rotating
		(int patch_overlap_width,
		 int ip_interp_handle, int ip_interp_par_table_handle,
		 bool print_msg_flag);
	void setup_ghost_zones__plus_xz_quadrant_mirrored
		(int patch_overlap_width,
		 int ip_interp_handle, int ip_interp_par_table_handle,
		 bool print_msg_flag);
	void setup_ghost_zones__plus_xz_quadrant_rotating
		(int patch_overlap_width,
		 int ip_interp_handle, int ip_interp_par_table_handle,
		 bool print_msg_flag);
	void setup_ghost_zones__plus_xyz_octant_mirrored
		(int patch_overlap_width,
		 int ip_interp_handle, int ip_interp_par_table_handle,
		 bool print_msg_flag);
	void setup_ghost_zones__plus_xyz_octant_rotating
		(int patch_overlap_width,
		 int ip_interp_handle, int ip_interp_par_table_handle,
		 bool print_msg_flag);

	// create/interlink a pair of periodic-symmetry ghost zones
	static
	  void create_periodic_symmetry_ghost_zones
		(const patch_edge& ex, const patch_edge& ey,
		 bool ipar_map_is_plus);

	// construct a pair of interpatch ghost zones
	// ... automagically figures out which edges are adjacent
	static
	  void create_interpatch_ghost_zones(patch &px, patch &py,
					     int patch_overlap_width);

	// finish setup of a pair of interpatch ghost zones
	// ... automagically figures out which edges are adjacent
	static
	  void finish_interpatch_setup
		(patch &px, patch &py,
		 int patch_overlap_width,
		 int ip_interp_handle, int ip_interp_par_table_handle);

	// assert() that all ghost zones of all patches are fully setup
	void assert_all_ghost_zones_fully_setup() const;


private:
        // we forbid copying and passing by value
        // by declaring the copy constructor and assignment operator
        // private, but never defining them
	patch_system(const patch_system &rhs);
	patch_system& operator=(const patch_system &rhs);

private:
	// local <--> global coordinate mapping
	global_coords global_coords_;

	// meta-info about patch system
	enum patch_system_type type_;
	int N_patches_;
	int N_grid_points_, ghosted_N_grid_points_;
        int max_N_additional_points_;
        int N_additional_points_;

	// [pn] = --> individual patches
	// *** constructor initialization list ordering:
	// *** this must be declared after  N_patches_
	patch** all_patches_;

	// [pn] = starting grid point number of individual patches
	// ... arrays are actually of size N_patches_+1, the [N_patches_]
	//     entries are == N_grid_points_ and ghosted_N_grid_points_
	// *** constructor initialization list ordering:
	// *** these must be declared after  N_patches_
	int* starting_gpn_;
	int* ghosted_starting_gpn_;

	// pointers to storage blocks for all gridfns
	// ... patches point into these, but we own the storage blocks
	fp* gridfn_storage_;
	fp* ghosted_gridfn_storage_;

	// min/max m over all ghost zone points
        mutable int global_min_ym_;
	mutable int global_max_ym_;

	// info about the surface interpolator
	// ... used only by radius_in_local_xyz_direction()
	int surface_interp_handle_, surface_interp_par_table_handle_;
	};

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif
