// patch_system.cc -- describes the (an) entire system of patches
// $Header$

//
// patch_system::patch_system
// patch_system::~patch_system
// patch_system::create_patches
// patch_system::setup_gridfn_storage
// patch_system::setup_ghost_zones__full_sphere
// patch_system::setup_ghost_zones__plus_z_hemisphere
// patch_system::setup_ghost_zones__plus_xy_quadrant_mirrored
// patch_system::setup_ghost_zones__plus_xy_quadrant_rotating
// patch_system::setup_ghost_zones__plus_xz_quadrant_mirrored
// patch_system::setup_ghost_zones__plus_xz_quadrant_rotating
// patch_system::setup_ghost_zones__plus_xyz_octant_mirrored
// patch_system::setup_ghost_zones__plus_xyz_octant_rotating
// patch_system::create_periodic_symmetry_ghost_zones
// patch_system::create_interpatch_ghost_zones
// patch_system::finish_interpatch_setup
// patch_system::assert_all_ghost_zones_fully_setup
//
// patch_system::N_patches_of_type
// patch_system::name_of_type
// patch_system::type_of_name
// patch_system::plus_or_minus_xyz_patch
// patch_system::patch_number_of_name
//
// patch_system::patch_irho_isigma_of_gpn
// patch_system::ghosted_patch_irho_isigma_of_gpn
//
// patch_system::set_gridfn_to_constant
// patch_system::gridfn_copy
// patch_system::ghosted_gridfn_copy
// patch_system::add_to_gridfn
// patch_system::add_to_ghosted_gridfn
// patch_system::gridfn_norms
// patch_system::ghosted_gridfn_norms
//
// patch_system::circumference
// patch_system::integrate_gridfn
// 
// patch_system::patch_containing_local_xyz
// patch_system::radius_in_local_xyz_direction
// patch_system::radii_in_local_xyz_directions
//
// patch_system::print_unknown_gridfn
// patch_system::read_unknown_gridfn
// patch_system::output_unknown_gridfn
//
// patch_system::synchronize
// patch_system::compute_synchronize_Jacobian
// patch_system::synchronize_Jacobian_global_minmax_ym
// patch_system::synchronize_Jacobian
/// patch_system::fold_Jacobian
/// patch_system::ghost_zone_Jacobian
//

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <vector>
#include <sstream>
#include <string>

#include "../cctk.h"



#include "../gr/gfns.hh"
#include "../jtutil/util.hh"
#include "../jtutil/array.hh"
#include "../jtutil/cpm_map.hh"
#include "../jtutil/linear_map.hh"

using namespace SAMRAI;

#ifdef HAVE_CAPABILITY_HDF5
// Some macros to fix compatibility issues as long as 1.8.0 is in beta
// phase
#  define H5_USE_16_API 1
#  include <hdf5.h>
#endif

#include "coords.hh"
#include "grid.hh"
#include "fd_grid.hh"
#include "patch.hh"
#include "patch_edge.hh"
#include "patch_interp.hh"
#include "ghost_zone.hh"
#include "patch_info.hh"
#include "patch_system.hh"
#include "patch_system_info.hh"

// all the code in this file is inside this namespace
namespace AHFD
	  {

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function constructs a  patch_system  object.
//
// Constructor arguments:
// ghost_zone_width = Width in grid points of all ghost zones.
// patch_overlap_width = Number of grid points that adjacent
//		      nominally-just-touching patches should overlap.
//		      For example, with patch_overlap_width == 3, here
//		      are the grid points of two neighboring patches:
//			x   x   x   x   x   X
//                                      |
//			            O   o   o   o   o   o
//		      Here | marks the "just touching" boundary,
//		      x and o the grid points before this extension,
//		      and X and O the extra grid points added by this
//		      extension.  For this example, the patch_extend_width
//		      parameter used by some other functions would
//		      be 1; in general
//			patch_overlap_width = 2*patch_extend_width + 1
// N_zones_per_right_angle = This sets the grid spacing (same in both
//			     directions) to 90.0 / N_zones_per_right_angle.
//			     It's a fatal error (error_exit()) if this
//			     doesn't evenly divide the grid sizes in both
//			     directions.
// ip_interp_handle = Cactus handle to the interpatch interpolation operator;
//		      this must be a 1-D interpolator
// ip_interp_par_table_handle = Cactus handle to the parameter table for the
//				interpatch interpolation operator
// surface_interp_handle = Cactus handle to the surface interpolation
//			   operator; this is optional, and is only used by
//			     radius_in_{local,global}_xyz_direction()
//			   If this is used, it must be a 2-D interpolator
// surface_interp_par_table_handle = Cactus handle to the parameter table
//				     for the surface interpolation operator;
//				     this is optional, and is only used by
//				       radius_in_local_xyz_direction()
// print_summary_msg_flag = true to print 2 lines of CCTK_VInfo messages
//				 giving the patch system type and origin
//			    false to skip this
// print_detailed_msg_flag = true to print extensive messages tracing the
//				  creation and initialization of various
//				  data structures
//			     false to skip this
//
patch_system::patch_system(fp origin_x_in, fp origin_y_in, fp origin_z_in,
			   enum patch_system_type type_in,
			   int ghost_zone_width, int patch_overlap_width,
			   int N_zones_per_right_angle,
			   int min_gfn_in, int max_gfn_in,
			   int ghosted_min_gfn_in, int ghosted_max_gfn_in,
			   int ip_interp_handle, int ip_interp_par_table_handle,
			   int surface_interp_handle_in,
			   int surface_interp_par_table_handle_in,
			   bool print_summary_msg_flag,
			   bool print_detailed_msg_flag)

	: global_coords_(origin_x_in, origin_y_in, origin_z_in),
	  type_(type_in),
	  N_patches_(N_patches_of_type(type_in)),
	  all_patches_(new patch*[N_patches_]),
	  starting_gpn_(new int[N_patches_+1]),
	  ghosted_starting_gpn_(new int[N_patches_+1]),
	  gridfn_storage_(NULL),		// set in setup_gridfn_storage()
	  ghosted_gridfn_storage_(NULL),	// set in setup_gridfn_storage()
	  global_min_ym_(0), global_max_ym_(0),
					// set in compute_synchronize_Jacobian()
	  surface_interp_handle_          (surface_interp_handle_in),
	  surface_interp_par_table_handle_(surface_interp_par_table_handle_in)
{
if (! jtutil::is_odd(patch_overlap_width))
   then error_exit(ERROR_EXIT,
"***** patch_system::patch_system(): implementation restriction:\n"
"        patch_overlap_width=%d, but we only support odd values!\n"
,
		   patch_overlap_width);			/*NOTREACHED*/
const int patch_extend_width = patch_overlap_width >> 1;

if (ghost_zone_width < fd_grid::molecule_radius())
   then error_exit(ERROR_EXIT,
"***** patch_system::patch_system():\n"
"        must have ghost_zone_width >= fd_grid::molecule_radius()\n"
"        but got ghost_zone_width=%d fd_grid::molecule_radius()=%d!\n"
"        FINITE_DIFF_ORDER=%d (see #define in src/include/config.hh)\n"
,
		   ghost_zone_width, fd_grid::molecule_radius(),
		   FINITE_DIFF_ORDER);				/*NOTREACHED*/

if (print_summary_msg_flag)
   then {
	CCTK_VInfo(CCTK_THORNSTRING,
		   "      constructing %s patch system",
		   name_of_type(type()));
	CCTK_VInfo(CCTK_THORNSTRING,
		   "                   at origin (%g,%g,%g)",
		   double(origin_x()), double(origin_y()), double(origin_z()));
	CCTK_VInfo(CCTK_THORNSTRING,
		   "                   with %d angular zones per right angle",
		   N_zones_per_right_angle);
	}
// construct/interlink the patches and ghost zones
switch	(type_in)
	{
case patch_system__full_sphere:
	create_patches(patch_system_info::full_sphere::patch_info_array,
		       ghost_zone_width, patch_extend_width,
		       N_zones_per_right_angle,
		       print_detailed_msg_flag);
	setup_gridfn_storage(min_gfn_in, max_gfn_in,
			     ghosted_min_gfn_in, ghosted_max_gfn_in,
			     print_detailed_msg_flag);
	setup_ghost_zones__full_sphere(patch_overlap_width,
				       ip_interp_handle,
				       ip_interp_par_table_handle,
				       print_detailed_msg_flag);
	break;

case patch_system__plus_z_hemisphere:
	create_patches(patch_system_info::plus_z_hemisphere::patch_info_array,
		       ghost_zone_width, patch_extend_width,
		       N_zones_per_right_angle,
		       print_detailed_msg_flag);
	setup_gridfn_storage(min_gfn_in, max_gfn_in,
			     ghosted_min_gfn_in, ghosted_max_gfn_in,
			     print_detailed_msg_flag);
	setup_ghost_zones__plus_z_hemisphere(patch_overlap_width,
					     ip_interp_handle,
					     ip_interp_par_table_handle,
					     print_detailed_msg_flag);
	break;

case patch_system__plus_xy_quadrant_mirrored:
	create_patches(patch_system_info::plus_xy_quadrant::patch_info_array,
		       ghost_zone_width, patch_extend_width,
		       N_zones_per_right_angle,
		       print_detailed_msg_flag);
	setup_gridfn_storage(min_gfn_in, max_gfn_in,
			     ghosted_min_gfn_in, ghosted_max_gfn_in,
			     print_detailed_msg_flag);
	setup_ghost_zones__plus_xy_quadrant_mirrored(patch_overlap_width,
						     ip_interp_handle,
						     ip_interp_par_table_handle,
						     print_detailed_msg_flag);
	break;

case patch_system__plus_xy_quadrant_rotating:
	create_patches(patch_system_info::plus_xy_quadrant::patch_info_array,
		       ghost_zone_width, patch_extend_width,
		       N_zones_per_right_angle,
		       print_detailed_msg_flag);
	setup_gridfn_storage(min_gfn_in, max_gfn_in,
			     ghosted_min_gfn_in, ghosted_max_gfn_in,
			     print_detailed_msg_flag);
	setup_ghost_zones__plus_xy_quadrant_rotating(patch_overlap_width,
						     ip_interp_handle,
						     ip_interp_par_table_handle,
						     print_detailed_msg_flag);
	break;

case patch_system__plus_xz_quadrant_mirrored:
	create_patches(patch_system_info::plus_xz_quadrant::patch_info_array,
		       ghost_zone_width, patch_extend_width,
		       N_zones_per_right_angle,
		       print_detailed_msg_flag);
	setup_gridfn_storage(min_gfn_in, max_gfn_in,
			     ghosted_min_gfn_in, ghosted_max_gfn_in,
			     print_detailed_msg_flag);
	setup_ghost_zones__plus_xz_quadrant_mirrored(patch_overlap_width,
						     ip_interp_handle,
						     ip_interp_par_table_handle,
						     print_detailed_msg_flag);
	break;

case patch_system__plus_xz_quadrant_rotating:
	create_patches(patch_system_info::plus_xz_quadrant::patch_info_array,
		       ghost_zone_width, patch_extend_width,
		       N_zones_per_right_angle,
		       print_detailed_msg_flag);
	setup_gridfn_storage(min_gfn_in, max_gfn_in,
			     ghosted_min_gfn_in, ghosted_max_gfn_in,
			     print_detailed_msg_flag);
	setup_ghost_zones__plus_xz_quadrant_rotating(patch_overlap_width,
						     ip_interp_handle,
						     ip_interp_par_table_handle,
						     print_detailed_msg_flag);
	break;

case patch_system__plus_xyz_octant_mirrored:
	create_patches(patch_system_info::plus_xyz_octant::patch_info_array,
		       ghost_zone_width, patch_extend_width,
		       N_zones_per_right_angle,
		       print_detailed_msg_flag);
	setup_gridfn_storage(min_gfn_in, max_gfn_in,
			     ghosted_min_gfn_in, ghosted_max_gfn_in,
			     print_detailed_msg_flag);
	setup_ghost_zones__plus_xyz_octant_mirrored(patch_overlap_width,
						    ip_interp_handle,
						    ip_interp_par_table_handle,
						    print_detailed_msg_flag);
	break;

case patch_system__plus_xyz_octant_rotating:
	create_patches(patch_system_info::plus_xyz_octant::patch_info_array,
		       ghost_zone_width, patch_extend_width,
		       N_zones_per_right_angle,
		       print_detailed_msg_flag);
	setup_gridfn_storage(min_gfn_in, max_gfn_in,
			     ghosted_min_gfn_in, ghosted_max_gfn_in,
			     print_detailed_msg_flag);
	setup_ghost_zones__plus_xyz_octant_rotating(patch_overlap_width,
						    ip_interp_handle,
						    ip_interp_par_table_handle,
						    print_detailed_msg_flag);
	break;

default:
	error_exit(ERROR_EXIT,
"***** patch_system::patch_system(): bad type_in=(int)%d!\n"
,
		   int(type_in));				/*NOTREACHED*/
	}

if (print_summary_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      ==> %d nominal, %d ghosted angular grid points",
		   N_grid_points(), ghosted_N_grid_points());
}

//******************************************************************************

//
// This function destroys a  patch_system  object.
//
patch_system::~patch_system()
{
delete[] ghosted_gridfn_storage_;
delete[] gridfn_storage_;
delete[] ghosted_starting_gpn_;
delete[] starting_gpn_;

	for (int pn = N_patches()-1 ; pn >= 0 ; --pn)
	{
	delete &ith_patch(pn);
	}

delete[] all_patches_;
}

//******************************************************************************

//
// This function is called from the patch_system:: constructor to
// construct a set of patches as described by an array of patch_info
// structures and associated arguments, and make these patches members
// of this patch system.  This function also correctly sets
//	N_grid_points_
//	N_ghosted_grid_points_
//	max_N_additional_points_
//	N_additional_points_
//	all_patches_[]
//	starting_gpn_[]
//	ghosted_starting_gpn_[]
// This function does *NOT* create any of the ghost zones, and does
// *NOT* set up any gridfns.
//
void patch_system::create_patches(const struct patch_info patch_info_in[],
				  int ghost_zone_width, int patch_extend_width,
				  int N_zones_per_right_angle,
				  bool print_msg_flag)
{
max_N_additional_points_ = 1;
// N_additional_points_ = 0;
N_additional_points_ = 1;
N_grid_points_ = 0;
ghosted_N_grid_points_ = 0;
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	const struct patch_info& pi = patch_info_in[pn];
	const struct grid::grid_array_pars& grid_array_pars
		= pi.grid_array_pars(ghost_zone_width,
				     patch_extend_width,
				     N_zones_per_right_angle);
	const struct grid::grid_pars& grid_pars
		= pi.grid_pars(patch_extend_width,
			       N_zones_per_right_angle);

	if (print_msg_flag)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "         constructing %s patch (%d x %d grid points)",
			   pi.name,
			   jtutil::how_many_in_range(grid_array_pars.min_irho,
						     grid_array_pars.max_irho),
			   jtutil::how_many_in_range(grid_array_pars.min_isigma,
						     grid_array_pars.max_isigma));

	struct patch *p;
	switch	(pi.ctype)
		{
	case 'z':
		p = new z_patch(*this, pn,
				pi.name, pi.is_plus,
				grid_array_pars, grid_pars);
		break;
	case 'x':
		p = new x_patch(*this, pn,
				pi.name, pi.is_plus,
				grid_array_pars, grid_pars);
		break;
	case 'y':
		p = new y_patch(*this, pn,
				pi.name, pi.is_plus,
				grid_array_pars, grid_pars);
		break;
	default:
		error_exit(ERROR_EXIT,
"***** patch_system::create_patches():\n"
"        unknown patch_info_in[pn=%d].ctype=0x%02d='%c'!\n"
,
			   pn, pi.ctype, pi.ctype);		/*NOTREACHED*/
		}

	// these record number of grid points in *previous* patches,
	// i.e. they do *not* include the number of grid points in this patch
	starting_gpn_[pn] = N_grid_points_;
	ghosted_starting_gpn_[pn] = ghosted_N_grid_points_;

	N_grid_points_+= p->N_grid_points();
	ghosted_N_grid_points_ += p->ghosted_N_grid_points();

	all_patches_[pn] = p;
	}

starting_gpn_        [N_patches_] = N_grid_points_;
ghosted_starting_gpn_[N_patches_] = ghosted_N_grid_points_;
}

//******************************************************************************

//
// This function is called from the patch_system:: constructor to set
// up the storage for all gridfns in all patches, giving each gridfn a
// contiguous-across-all-patches storage array.  This function also makes
// a number of self-consistency checks to ensure that the gridfn storage
// subscripting is set up properly.
//
// This function assumes that all the patches have already been constructed
// before it is called.
//
// For example, given the patches {x,y,z}, the ghosted gridfns {H,J},
// and the nominal gridfns {a,b,c}, we might picture the storage like
// this:
//
//	xa xa xa ya ya za za za za
//	xb xb xb yb yb zb zb zb zb
//	xc xc xc yc yc zc zc zc zc
//
//	xH xH xH xH yH yH yH zH zH zH zH zH
//	xJ xJ xJ xJ yJ yJ yJ zJ zJ zJ zJ zJ
//
// Here the upper/lower blocks are for nominal/ghosted  gridfns.
// The storage is taken as being contiguous within each row (in fact
// within each block).  Thus the storage for all the nominal gridfns
// (or all the ghosted gridfns) in a single patch is *non*-contiguous.
//
// The creation of patches is done in several phases: first the patches
// are constructed with no gridfn storage, then we are called to set up
// the gridfn storage (taking into account the sizes of the other patches),
// then finally ghost zones are constructed and interlinked.
//
// FIXME: We should pad the gridfn storage as necessary to avoid cache
//        conflicts, but we don't do this at present.
//
void patch_system::setup_gridfn_storage
	(int min_gfn_in, int max_gfn_in,
	 int ghosted_min_gfn_in, int ghosted_max_gfn_in,
	 bool print_msg_flag)
{
const int N_gridfns_in = jtutil::how_many_in_range(min_gfn_in, max_gfn_in);
const int ghosted_N_gridfns_in
	= jtutil::how_many_in_range(ghosted_min_gfn_in, ghosted_max_gfn_in);

const int gfn_stride = N_grid_points() + max_N_additional_points();
const int ghosted_gfn_stride = ghosted_N_grid_points() + max_N_additional_points();

const int N_storage = gfn_stride * N_gridfns_in;
const int ghosted_N_storage = ghosted_gfn_stride * ghosted_N_gridfns_in;

if (print_msg_flag)
   then {
	CCTK_VInfo(CCTK_THORNSTRING,
		   "      setting up gridfn storage");
	CCTK_VInfo(CCTK_THORNSTRING,
		   "         gfn=[%d,%d] ghosted_gfn=[%d,%d]",
		   min_gfn_in, max_gfn_in,
		   ghosted_min_gfn_in, ghosted_max_gfn_in);
	CCTK_VInfo(CCTK_THORNSTRING,
		   "         N_grid_points()=%d ghosted_N_grid_points()=%d",
		   N_grid_points(), ghosted_N_grid_points());
	}

// storage arrays for all gridfns
gridfn_storage_         = new fp[N_storage];
ghosted_gridfn_storage_ = new fp[ghosted_N_storage];

// divide up the storage array among the patches
// and set up the storage in the individual patches themselves
	  {
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	const int         posn =         starting_gpn_[pn];
	const int ghosted_posn = ghosted_starting_gpn_[pn];
	const struct grid_arrays::gridfn_pars gridfn_pars
		= {
		  min_gfn_in, max_gfn_in,
		  & gridfn_storage_[posn],
		  gfn_stride, 0, 0
		  };
	const struct grid_arrays::gridfn_pars ghosted_gridfn_pars
		= {
		  ghosted_min_gfn_in, ghosted_max_gfn_in,
		  & ghosted_gridfn_storage_[ghosted_posn],
		  ghosted_gfn_stride, 0, 0
		  };

	patch& p = ith_patch(pn);
	p.setup_gridfn_storage(gridfn_pars, ghosted_gridfn_pars);
	}
	  }

if (print_msg_flag)
   then {
	CCTK_VInfo(CCTK_THORNSTRING,
		   "         checking that storage is partitioned properly");
	}

// check to make sure storage for distinct gridfns
// forms a partition of the overall storage array
const patch& pfirst = ith_patch(0);
const patch& plast  = ith_patch(N_patches()-1);
	  {
	for (int gfn = min_gfn() ; gfn+1 < max_gfn() ; ++gfn)
	{
	// range of storage occupied by gridfns:
	//	gfn   --> [gfn_first, gfn_last]
	//	gfn+1 --> [gfn1_first, gfn1_last]
	const fp* const gfn_last_ptr
		= & plast.gridfn(gfn, plast.max_irho(),
				      plast.max_isigma());
	const fp* const gfn1_first_ptr
		= & pfirst.gridfn(gfn+1, pfirst.min_irho(),
					 pfirst.min_isigma());
	if (! (gfn1_first_ptr == gfn_last_ptr+max_N_additional_points()+1))
	   then error_exit(PANIC_EXIT,
"***** patch_system::setup_gridfn_storage():\n"
"        nominal-grid gridfns don't partition overall storage array!"
"        (this should never happen!)\n"
"        gfn=%d last point at gfn_last_ptr=%p\n"
"	 gfn+1=%d first point at gfn1_first_ptr=%p\n"
"        should have gfn1_first_ptr == gfn_last_ptr+1\n"
,
			   gfn, (const void *) gfn_last_ptr,
			   gfn+1, (const void *) gfn1_first_ptr); /*NOTREACHED*/
	}
	  }

	  {
	for (int ghosted_gfn = ghosted_min_gfn() ;
	     ghosted_gfn+1 < ghosted_max_gfn() ;
	     ++ghosted_gfn)
	{
	// range of storage occupied by ghosted gridfns:
	//	ghosted_gfn   --> [gfn_first, gfn_last]
	//	ghosted_gfn+1 --> [gfn1_first, gfn1_last]
	const fp* const ghosted_gfn_last_ptr
		= & plast.ghosted_gridfn(ghosted_gfn,
					 plast.ghosted_max_irho(),
					 plast.ghosted_max_isigma());
	const fp* const ghosted_gfn1_first_ptr
		= & pfirst.ghosted_gridfn(ghosted_gfn+1,
					  pfirst.ghosted_min_irho(),
					  pfirst.ghosted_min_isigma());
	if (! (ghosted_gfn1_first_ptr == ghosted_gfn_last_ptr+max_N_additional_points()+1))
	   then error_exit(PANIC_EXIT,
"***** patch_system::setup_gridfn_storage():\n"
"        ghosted-grid gridfns don't partition overall storage array!"
"        (this should never happen!)\n"
"        ghosted_gfn=%d last point at ghosted_gfn_last_ptr=%p\n"
"	 ghosted_gfn+1=%d first point at ghosted_gfn1_first_ptr=%p\n"
"        should have ghosted_gfn1_first_ptr == ghosted_gfn_last_ptr+1\n"
,
			   ghosted_gfn, (const void *) ghosted_gfn_last_ptr,
			   ghosted_gfn+1,
			      (const void *) ghosted_gfn1_first_ptr);
								/*NOTREACHED*/
	}
	  }

// check to make sure storage for distinct patches
// forms a partition of the storage for each gridfn
	  {
	for (int gfn = min_gfn() ; gfn < max_gfn() ; ++gfn)
	{
		for (int pn = 0 ; pn+1 < N_patches() ; ++pn)
		{
		const patch& p = ith_patch(pn);
		const patch& p1 = ith_patch(pn+1);

		// range of storage occupied by gridfn:
		//	p  --> [p_first, p_last]
		//	p1 --> [p1_first, p1_last]
		const fp* const p_last_ptr
			= & p.gridfn(gfn, p.max_irho(), p.max_isigma());
		const fp* const p1_first_ptr
			= & p1.gridfn(gfn, p1.min_irho(), p1.min_isigma());
		if (! (p1_first_ptr == p_last_ptr+1))
		   then error_exit(PANIC_EXIT,
"***** patch_system::setup_gridfn_storage():\n"
"        nominal-grid patches gridfns don't partition storage for gfn=%d!\n"
"        (this should never happen!)\n"
"        gfn=%d %s patch last point at p_last_ptr=%p\n"
"	 gfn=%d %s patch first point at p1_first_ptr=%p\n"
"        should have p1_first_ptr == p_last_ptr+1\n"
,
			   gfn,
			   gfn, p.name(), (const void *) p_last_ptr,
			   gfn+1, p1.name(), (const void *) p1_first_ptr);
								 /*NOTREACHED*/
		}
	}
	  }

	  {
	for (int ghosted_gfn = ghosted_min_gfn() ;
	     ghosted_gfn < ghosted_max_gfn() ;
	     ++ghosted_gfn)
	{
		for (int pn = 0 ; pn+1 < N_patches() ; ++pn)
		{
		const patch& p = ith_patch(pn);
		const patch& p1 = ith_patch(pn+1);

		// range of storage occupied by ghosted gridfn:
		//	p  --> [p_first, p_last]
		//	p1 --> [p1_first, p1_last]
		const fp* const p_last_ptr
			= & p.ghosted_gridfn(ghosted_gfn,
					     p.ghosted_max_irho(),
					     p.ghosted_max_isigma());
		const fp* const p1_first_ptr
			= & p1.ghosted_gridfn(ghosted_gfn,
					      p1.ghosted_min_irho(),
					      p1.ghosted_min_isigma());
		if (! (p1_first_ptr == p_last_ptr+1))
		   then error_exit(PANIC_EXIT,
"***** patch_system::setup_gridfn_storage():\n"
"        ghosted-grid patches gridfns don't partition storage for gfn=%d!\n"
"        (this should never happen!)\n"
"        %s patch (pn=%d) last point at p_last_ptr=%p\n"
"	 %s patch (pn=%d) first point at p1_first_ptr=%p\n"
"        should have p1_first_ptr == p_last_ptr+1\n"
,
			   ghosted_gfn,
			   p.name(), pn, (const void *) p_last_ptr,
			   p1.name(), pn+1, (const void *) p1_first_ptr);
								 /*NOTREACHED*/
		}
	}
	  }
}

//******************************************************************************

//
// This function sets up (constructs and interlinks) the ghost zones
// for a full-sphere patch system.
//
void patch_system::setup_ghost_zones__full_sphere
	(int patch_overlap_width,
	 int ip_interp_handle, int ip_interp_par_table_handle,
	 bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      seting up full sphere ghost zones");

patch& pz = ith_patch(patch_system_info::full_sphere::patch_number__pz);
patch& px = ith_patch(patch_system_info::full_sphere::patch_number__px);
patch& py = ith_patch(patch_system_info::full_sphere::patch_number__py);
patch& mx = ith_patch(patch_system_info::full_sphere::patch_number__mx);
patch& my = ith_patch(patch_system_info::full_sphere::patch_number__my);
patch& mz = ith_patch(patch_system_info::full_sphere::patch_number__mz);
// create the ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         creating ghost zones");
create_interpatch_ghost_zones(pz, px, patch_overlap_width);
create_interpatch_ghost_zones(pz, py, patch_overlap_width);
create_interpatch_ghost_zones(pz, mx, patch_overlap_width);
create_interpatch_ghost_zones(pz, my, patch_overlap_width);
create_interpatch_ghost_zones(px, py, patch_overlap_width);
create_interpatch_ghost_zones(py, mx, patch_overlap_width);
create_interpatch_ghost_zones(mx, my, patch_overlap_width);
create_interpatch_ghost_zones(my, px, patch_overlap_width);
create_interpatch_ghost_zones(mz, px, patch_overlap_width);
create_interpatch_ghost_zones(mz, py, patch_overlap_width);
create_interpatch_ghost_zones(mz, mx, patch_overlap_width);
create_interpatch_ghost_zones(mz, my, patch_overlap_width);

// finish setting up the interpatch ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         finishing interpatch setup");
finish_interpatch_setup(pz, px,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(pz, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(pz, mx,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(pz, my,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(px, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(py, mx,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(mx, my,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(my, px,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
// MARKS CRASH
finish_interpatch_setup(mz, px,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(mz, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(mz, mx,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(mz, my,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);

assert_all_ghost_zones_fully_setup();
}

//******************************************************************************

//
// This function sets up (constructs and interlinks) the ghost zones
// for a +z hemisphere patch system.
//
void patch_system::setup_ghost_zones__plus_z_hemisphere
	(int patch_overlap_width,
	 int ip_interp_handle, int ip_interp_par_table_handle,
	 bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      setting up +z hemisphere ghost zones");

patch& pz = ith_patch(patch_system_info::plus_z_hemisphere::patch_number__pz);
patch& px = ith_patch(patch_system_info::plus_z_hemisphere::patch_number__px);
patch& py = ith_patch(patch_system_info::plus_z_hemisphere::patch_number__py);
patch& mx = ith_patch(patch_system_info::plus_z_hemisphere::patch_number__mx);
patch& my = ith_patch(patch_system_info::plus_z_hemisphere::patch_number__my);

// create the ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         creating ghost zones");
create_interpatch_ghost_zones(pz, px, patch_overlap_width);
create_interpatch_ghost_zones(pz, py, patch_overlap_width);
create_interpatch_ghost_zones(pz, mx, patch_overlap_width);
create_interpatch_ghost_zones(pz, my, patch_overlap_width);
create_interpatch_ghost_zones(px, py, patch_overlap_width);
create_interpatch_ghost_zones(py, mx, patch_overlap_width);
create_interpatch_ghost_zones(mx, my, patch_overlap_width);
create_interpatch_ghost_zones(my, px, patch_overlap_width);
px.create_mirror_symmetry_ghost_zone(px.max_rho_patch_edge());
py.create_mirror_symmetry_ghost_zone(py.max_rho_patch_edge());
mx.create_mirror_symmetry_ghost_zone(mx.min_rho_patch_edge());
my.create_mirror_symmetry_ghost_zone(my.min_rho_patch_edge());

// finish setting up the interpatch ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         finishing interpatch setup");
finish_interpatch_setup(pz, px,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(pz, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(pz, mx,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(pz, my,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(px, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(py, mx,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(mx, my,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(my, px,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);

assert_all_ghost_zones_fully_setup();
}

//******************************************************************************

//
// This function sets up (constructs and interlinks) the ghost zones
// for a +xy quadrant (mirrored) patch system.
//
void patch_system::setup_ghost_zones__plus_xy_quadrant_mirrored
	(int patch_overlap_width,
	 int ip_interp_handle, int ip_interp_par_table_handle,
	 bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      setting up +xy quadrant (mirrored) ghost zones");

patch& pz = ith_patch(patch_system_info::plus_xy_quadrant::patch_number__pz);
patch& px = ith_patch(patch_system_info::plus_xy_quadrant::patch_number__px);
patch& py = ith_patch(patch_system_info::plus_xy_quadrant::patch_number__py);
patch& mz = ith_patch(patch_system_info::plus_xy_quadrant::patch_number__mz);

// create the ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         creating ghost zones");
create_interpatch_ghost_zones(pz, px, patch_overlap_width);
create_interpatch_ghost_zones(pz, py, patch_overlap_width);
create_interpatch_ghost_zones(px, py, patch_overlap_width);
create_interpatch_ghost_zones(mz, px, patch_overlap_width);
create_interpatch_ghost_zones(mz, py, patch_overlap_width);
pz.create_mirror_symmetry_ghost_zone(pz.min_rho_patch_edge());
pz.create_mirror_symmetry_ghost_zone(pz.min_sigma_patch_edge());
px.create_mirror_symmetry_ghost_zone(px.min_sigma_patch_edge());
py.create_mirror_symmetry_ghost_zone(py.max_sigma_patch_edge());
mz.create_mirror_symmetry_ghost_zone(mz.max_rho_patch_edge());
mz.create_mirror_symmetry_ghost_zone(mz.max_sigma_patch_edge());

// finish setting up the interpatch ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         finishing interpatch setup");
finish_interpatch_setup(pz, px,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(pz, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(px, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(mz, px,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(mz, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);

assert_all_ghost_zones_fully_setup();
}

//******************************************************************************

//
// This function sets up (constructs and interlinks) the ghost zones
// for a +xy quadrant (rotating) patch system.
//
void patch_system::setup_ghost_zones__plus_xy_quadrant_rotating
	(int patch_overlap_width,
	 int ip_interp_handle, int ip_interp_par_table_handle,
	 bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      setting up +xy quadrant (rotating) ghost zones");

patch& pz = ith_patch(patch_system_info::plus_xy_quadrant::patch_number__pz);
patch& px = ith_patch(patch_system_info::plus_xy_quadrant::patch_number__px);
patch& py = ith_patch(patch_system_info::plus_xy_quadrant::patch_number__py);
patch& mz = ith_patch(patch_system_info::plus_xy_quadrant::patch_number__mz);

// create the ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         creating ghost zones");
create_interpatch_ghost_zones(pz, px, patch_overlap_width);
create_interpatch_ghost_zones(pz, py, patch_overlap_width);
create_interpatch_ghost_zones(px, py, patch_overlap_width);
create_interpatch_ghost_zones(mz, px, patch_overlap_width);
create_interpatch_ghost_zones(mz, py, patch_overlap_width);
create_periodic_symmetry_ghost_zones(pz.min_rho_patch_edge(),
				     pz.min_sigma_patch_edge(),
				     true);
create_periodic_symmetry_ghost_zones(px.min_sigma_patch_edge(),
				     py.max_sigma_patch_edge(),
				     true);
create_periodic_symmetry_ghost_zones(mz.max_rho_patch_edge(),
				     mz.max_sigma_patch_edge(),
				     true);

// finish setting up the interpatch ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         finishing interpatch setup");
finish_interpatch_setup(pz, px,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(pz, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(px, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(mz, px,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(mz, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);

assert_all_ghost_zones_fully_setup();
}

//******************************************************************************

//
// This function sets up (constructs and interlinks) the ghost zones
// for a +xz quadrant (mirrored) patch system.
//
void patch_system::setup_ghost_zones__plus_xz_quadrant_mirrored
	(int patch_overlap_width,
	 int ip_interp_handle, int ip_interp_par_table_handle,
	 bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      setting up +xz quadrant (mirrored) ghost zones");

patch& pz = ith_patch(patch_system_info::plus_xz_quadrant::patch_number__pz);
patch& px = ith_patch(patch_system_info::plus_xz_quadrant::patch_number__px);
patch& py = ith_patch(patch_system_info::plus_xz_quadrant::patch_number__py);
patch& my = ith_patch(patch_system_info::plus_xz_quadrant::patch_number__my);

// create the ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         creating ghost zones");
create_interpatch_ghost_zones(pz, px, patch_overlap_width);
create_interpatch_ghost_zones(pz, py, patch_overlap_width);
create_interpatch_ghost_zones(pz, my, patch_overlap_width);
create_interpatch_ghost_zones(px, py, patch_overlap_width);
create_interpatch_ghost_zones(px, my, patch_overlap_width);
pz.create_mirror_symmetry_ghost_zone(pz.min_sigma_patch_edge());
px.create_mirror_symmetry_ghost_zone(px.max_rho_patch_edge());
py.create_mirror_symmetry_ghost_zone(py.max_rho_patch_edge());
py.create_mirror_symmetry_ghost_zone(py.max_sigma_patch_edge());
my.create_mirror_symmetry_ghost_zone(my.min_rho_patch_edge());
my.create_mirror_symmetry_ghost_zone(my.min_sigma_patch_edge());

// finish setting up the interpatch ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         finishing interpatch setup");
finish_interpatch_setup(pz, px,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(pz, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(pz, my,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(px, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(px, my,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);

assert_all_ghost_zones_fully_setup();
}

//******************************************************************************

//
// This function sets up (constructs and interlinks) the ghost zones
// for a +xz quadrant (rotating) patch system.
//
void patch_system::setup_ghost_zones__plus_xz_quadrant_rotating
	(int patch_overlap_width,
	 int ip_interp_handle, int ip_interp_par_table_handle,
	 bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      setting up +xz quadrant (rotating) ghost zones");

patch& pz = ith_patch(patch_system_info::plus_xz_quadrant::patch_number__pz);
patch& px = ith_patch(patch_system_info::plus_xz_quadrant::patch_number__px);
patch& py = ith_patch(patch_system_info::plus_xz_quadrant::patch_number__py);
patch& my = ith_patch(patch_system_info::plus_xz_quadrant::patch_number__my);

// create the ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         creating ghost zones");
create_interpatch_ghost_zones(pz, px, patch_overlap_width);
create_interpatch_ghost_zones(pz, py, patch_overlap_width);
create_interpatch_ghost_zones(pz, my, patch_overlap_width);
create_interpatch_ghost_zones(px, py, patch_overlap_width);
create_interpatch_ghost_zones(px, my, patch_overlap_width);
px.create_mirror_symmetry_ghost_zone(px.max_rho_patch_edge());
py.create_mirror_symmetry_ghost_zone(py.max_rho_patch_edge());
my.create_mirror_symmetry_ghost_zone(my.min_rho_patch_edge());
create_periodic_symmetry_ghost_zones(pz.min_sigma_patch_edge(),
				     pz.min_sigma_patch_edge(),
				     false);
create_periodic_symmetry_ghost_zones(py.max_sigma_patch_edge(),
				     my.min_sigma_patch_edge(),
				     false);

// finish setting up the interpatch ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         finishing interpatch setup");
finish_interpatch_setup(pz, px,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(pz, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(pz, my,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(px, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(px, my,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);

assert_all_ghost_zones_fully_setup();
}

//******************************************************************************

//
// This function sets up (constructs and interlinks) the ghost zones
// for a +xyz octant (mirrored) patch system.
//
void patch_system::setup_ghost_zones__plus_xyz_octant_mirrored
	(int patch_overlap_width,
	 int ip_interp_handle, int ip_interp_par_table_handle,
	 bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      setting up +xyz octant (mirrored) ghost zones");

patch& pz = ith_patch(patch_system_info::plus_xyz_octant::patch_number__pz);
patch& px = ith_patch(patch_system_info::plus_xyz_octant::patch_number__px);
patch& py = ith_patch(patch_system_info::plus_xyz_octant::patch_number__py);

// create the ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         creating ghost zones");
create_interpatch_ghost_zones(pz, px, patch_overlap_width);
create_interpatch_ghost_zones(pz, py, patch_overlap_width);
create_interpatch_ghost_zones(px, py, patch_overlap_width);
pz.create_mirror_symmetry_ghost_zone(pz.min_rho_patch_edge());
pz.create_mirror_symmetry_ghost_zone(pz.min_sigma_patch_edge());
px.create_mirror_symmetry_ghost_zone(px.max_rho_patch_edge());
px.create_mirror_symmetry_ghost_zone(px.min_sigma_patch_edge());
py.create_mirror_symmetry_ghost_zone(py.max_rho_patch_edge());
py.create_mirror_symmetry_ghost_zone(py.max_sigma_patch_edge());

// finish setting up the interpatch ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         finishing interpatch setup");
finish_interpatch_setup(pz, px,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(pz, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(px, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);

assert_all_ghost_zones_fully_setup();
}

//******************************************************************************

//
// This function sets up (constructs and interlinks) the ghost zones
// for a +xyz octant (rotating) patch system.
//
void patch_system::setup_ghost_zones__plus_xyz_octant_rotating
	(int patch_overlap_width,
	 int ip_interp_handle, int ip_interp_par_table_handle,
	 bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      setting up +xyz octant (rotating) ghost zones");

patch& pz = ith_patch(patch_system_info::plus_xyz_octant::patch_number__pz);
patch& px = ith_patch(patch_system_info::plus_xyz_octant::patch_number__px);
patch& py = ith_patch(patch_system_info::plus_xyz_octant::patch_number__py);

// create the ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         creating ghost zones");
create_interpatch_ghost_zones(pz, px, patch_overlap_width);
create_interpatch_ghost_zones(pz, py, patch_overlap_width);
create_interpatch_ghost_zones(px, py, patch_overlap_width);
px.create_mirror_symmetry_ghost_zone(px.max_rho_patch_edge());
py.create_mirror_symmetry_ghost_zone(py.max_rho_patch_edge());
create_periodic_symmetry_ghost_zones(pz.min_rho_patch_edge(),
				     pz.min_sigma_patch_edge(),
				     true);
create_periodic_symmetry_ghost_zones(px.min_sigma_patch_edge(),
				     py.max_sigma_patch_edge(),
				     true);

// finish setting up the interpatch ghost zones
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         finishing interpatch setup");
finish_interpatch_setup(pz, px,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(pz, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);
finish_interpatch_setup(px, py,
			patch_overlap_width,
			ip_interp_handle, ip_interp_par_table_handle);

assert_all_ghost_zones_fully_setup();
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function creates a pair of periodic-symmetry ghost zones.
//
//static
  void patch_system::create_periodic_symmetry_ghost_zones
	(const patch_edge& ex, const patch_edge& ey,
	 bool ipar_map_is_plus)
{
ex.my_patch()
  .create_periodic_symmetry_ghost_zone(ex, ey, ipar_map_is_plus);

if (ex == ey)
   then {
	// ex and ey are the same edge (i.e. the symmetry maps the edge
	// back to itself), so we only want to set up the edge once
	// ==> no-op here
	}
   else ey.my_patch()
	  .create_periodic_symmetry_ghost_zone(ey, ex, ipar_map_is_plus);
}

//******************************************************************************

//
// This function automagically figures out which edges of two adjacent
// patches are adjacent, then creates both patches' ghost zones on those
// edges and interlinks them with their respective patches.
//
//static
  void patch_system::create_interpatch_ghost_zones(patch &px, patch &py,
						   int patch_overlap_width)
{
const patch_edge& ex = px.edge_adjacent_to_patch(py, patch_overlap_width);
const patch_edge& ey = py.edge_adjacent_to_patch(px, patch_overlap_width);

px.create_interpatch_ghost_zone(ex, ey, patch_overlap_width);
py.create_interpatch_ghost_zone(ey, ex, patch_overlap_width);
}

//******************************************************************************

//
// This function automagically figures out which edges of two adjacent
// patches are adjacent, then finishes setting up both ghost zones.
//
//static
  void patch_system::finish_interpatch_setup
	(patch &px, patch &py,
	 int patch_overlap_width,
	 int ip_interp_handle, int ip_interp_par_table_handle)
{
const patch_edge& ex = px.edge_adjacent_to_patch(py, patch_overlap_width);
const patch_edge& ey = py.edge_adjacent_to_patch(px, patch_overlap_width);

px.ghost_zone_on_edge(ex)
  .cast_to_interpatch_ghost_zone()
  .finish_setup(ip_interp_handle, ip_interp_par_table_handle);
py.ghost_zone_on_edge(ey)
  .cast_to_interpatch_ghost_zone()
  .finish_setup(ip_interp_handle, ip_interp_par_table_handle);
}

//******************************************************************************

//
// This function assert()s that all ghost zones of all patches have
// been fully set up.
//
void patch_system::assert_all_ghost_zones_fully_setup() const
{
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	ith_patch(pn).assert_all_ghost_zones_fully_setup();
	}
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function decodes a patch system's type into N_patches.
//
//static
  int patch_system::N_patches_of_type(enum patch_system_type type_in)
{
switch	(type_in)
	{
case patch_system__full_sphere:
	return patch_system_info::full_sphere::N_patches;
case patch_system__plus_z_hemisphere:
	return patch_system_info::plus_z_hemisphere::N_patches;
case patch_system__plus_xy_quadrant_mirrored:
case patch_system__plus_xy_quadrant_rotating:
	return patch_system_info::plus_xy_quadrant::N_patches;
case patch_system__plus_xz_quadrant_mirrored:
case patch_system__plus_xz_quadrant_rotating:
	return patch_system_info::plus_xz_quadrant::N_patches;
case patch_system__plus_xyz_octant_mirrored:
case patch_system__plus_xyz_octant_rotating:
	return patch_system_info::plus_xyz_octant::N_patches;
default:
	error_exit(PANIC_EXIT,
"***** patch_system::N_patches_of_type(): bad type=(int)%d!\n"
"                                         (this should never happen!)\n",
		   int(type_in));				/*NOTREACHED*/
	}
}

//******************************************************************************

//
// This function decodes a patch system's type into a human-readable
// type name (a C string).
//
//static
  const char* patch_system::name_of_type(enum patch_system_type type_in)
{
switch	(type_in)
	{
case patch_system__full_sphere:               return "full sphere";
case patch_system__plus_z_hemisphere:         return "+z hemisphere";
case patch_system__plus_xy_quadrant_mirrored: return "+xy quadrant (mirrored)";
case patch_system__plus_xy_quadrant_rotating: return "+xy quadrant (rotating)";
case patch_system__plus_xz_quadrant_mirrored: return "+xz quadrant (mirrored)";
case patch_system__plus_xz_quadrant_rotating: return "+xz quadrant (rotating)";
case patch_system__plus_xyz_octant_mirrored:  return "+xyz octant (mirrored)";
case patch_system__plus_xyz_octant_rotating:  return "+xyz octant (rotating)";
default:
	error_exit(PANIC_EXIT,
"***** patch_system::name_of_type(): bad type=(int)%d!\n"
"                                    (this should never happen!)\n",
		   int(type_in));				/*NOTREACHED*/
	}
}

//******************************************************************************

//
// This function encodes a human-readable type name (a C string) into
// a patch system's type into.
//
//static
  enum patch_system::patch_system_type
    patch_system::type_of_name(const char* name_in)
{
if	(STRING_EQUAL(name_in, "full sphere"))
   then return patch_system__full_sphere;
else if (STRING_EQUAL(name_in, "+z hemisphere"))
   then return patch_system__plus_z_hemisphere;
else if (STRING_EQUAL(name_in, "+xy quadrant (mirrored)"))
   then return patch_system__plus_xy_quadrant_mirrored;
else if (STRING_EQUAL(name_in, "+xy quadrant (rotating)"))
   then return patch_system__plus_xy_quadrant_rotating;
else if (STRING_EQUAL(name_in, "+xz quadrant (mirrored)"))
   then return patch_system__plus_xz_quadrant_mirrored;
else if (STRING_EQUAL(name_in, "+xz quadrant (rotating)"))
   then return patch_system__plus_xz_quadrant_rotating;
else if (STRING_EQUAL(name_in, "+xyz octant (mirrored)"))
   then return patch_system__plus_xyz_octant_mirrored;
else if (STRING_EQUAL(name_in, "+xyz octant (rotating)"))
   then return patch_system__plus_xyz_octant_rotating;
else	error_exit(PANIC_EXIT,
"***** patch_system::type_of_name(): unknown name=\"%s\"!",
		   name_in);					/*NOTREACHED*/
}

//******************************************************************************

//
// This function finds a (the) patch with a specified sign and xyz ctype.
// If no such patch exists, it does an  error_exit()  (and doesn't return
// to the caller).
//
// FIXME:
// - This function could be implemented to be very fast (using the
//   patch numbers in patch_system_info::), but right now it just does
//   a sequential search through all the patches, so it's pretty slow :(
//
const patch& patch_system::plus_or_minus_xyz_patch(bool is_plus, char ctype)
	const
{
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	const patch& p = ith_patch(pn);
	if ((p.is_plus() == is_plus) && (p.ctype() == ctype))
	   then return p;
	}

error_exit(ERROR_EXIT,
"***** patch_system::plus_or_minus_xyz_patch():\n"
"        can't find any %c%c patch!"
,
	   (is_plus ? '+' : '-'), ctype);			/*NOTREACHED*/
}

//******************************************************************************

//
// This function finds a patch from its human-readable name, and returns
// the patch number, or does an error_exit() if no patch is found with
// the specified name.
//
int patch_system::patch_number_of_name(const char* name) const
{
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	if (STRING_EQUAL(ith_patch(pn).name(), name))
	   then return pn;
	}

error_exit(ERROR_EXIT,
"***** patch_system::patch_number_of_name():\n"
"        no patch with name \"%s\"!\n",
	   name);						/*NOTREACHED*/
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function decodes a 0-origin grid point number into a
// (patch,irho,isigma) triple.
//
// Arguments:
// gpn = The grid point number to decode.
// (irho,isigma) = (out) The decoded patch coordinates.
//
// Results:
// This function returns a reference to the decoded patch.  (An alternative
// design would be to return this via a  patch*&  argument, but the design
// used here seems slightly cleaner to use in practice.)
//
const patch&
  patch_system::patch_irho_isigma_of_gpn(int gpn, int& irho, int& isigma)
	const
{
assert( gpn >= 0 );
assert( gpn < N_grid_points() );

	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	// n.b. [pn+1] is ok since starting_gpn_[] has size N_patches()+1
	if ((gpn >= starting_gpn_[pn])  && (gpn < starting_gpn_[pn+1]))
	   then {
		const patch& p = ith_patch(pn);
		const int gpn_in_patch = gpn - starting_gpn_[pn];
		assert( gpn_in_patch >= 0 );
		assert( gpn_in_patch < p.N_grid_points() );
		p.irho_isigma_of_gpn(gpn_in_patch, irho,isigma);
		return p;
		}
	}

error_exit(PANIC_EXIT,
"***** patch_system::patch_irho_isigma_of_gpn(gpn=%d):\n"
"        couldn't find any patch! (this should never happen!)\n"
"        N_grid_points()=%d\n"
,
	   gpn,
	   N_grid_points());					/*NOTREACHED*/
}

//******************************************************************************

//
// This function decodes a 0-origin grid point number into a *ghosted*
// (patch,irho,isigma) triple.
//
// Arguments:
// gpn = The grid point number to decode.
// (irho,isigma) = (out) The decoded patch coordinates.
//
// Results:
// This function returns a reference to the decoded patch.  (An alternative
// design would be to return this via a  patch*&  argument, but the design
// used here seems slightly cleaner to use in practice.)
//
const patch&
  patch_system::ghosted_patch_irho_isigma_of_gpn
	(int gpn, int& irho, int& isigma)
	const
{
assert( gpn >= 0 );
assert( gpn < ghosted_N_grid_points() );

	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	// n.b. [pn+1] is ok since ghosted_starting_gpn_[]
	//      has size N_patches()+1
	if (    (gpn >= ghosted_starting_gpn_[pn])
	     && (gpn <  ghosted_starting_gpn_[pn+1]))
	   then {
		const patch& p = ith_patch(pn);
		const int gpn_in_patch = gpn - ghosted_starting_gpn_[pn];
		assert( gpn_in_patch >= 0 );
		assert( gpn_in_patch < p.ghosted_N_grid_points() );
		p.ghosted_irho_isigma_of_gpn(gpn_in_patch, irho,isigma);
		return p;
		}
	}

error_exit(PANIC_EXIT,
"***** patch_system::ghosted_patch_irho_isigma_of_gpn(gpn=%d):\n"
"        couldn't find any patch! (this should never happen!)\n"
"        ghosted_N_grid_points()=%d\n"
,
	   gpn,
	   ghosted_N_grid_points());				/*NOTREACHED*/
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function sets a (nominal-grid) gridfn to a constant value.
//
void patch_system::set_gridfn_to_constant(fp a, int dst_gfn)
{
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	patch& p = ith_patch(pn);
		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		p.gridfn(dst_gfn, irho,isigma) = a;
		}
		}
	}
}

//******************************************************************************

//
// This function sets a ghosted gridfn to a constant value.
//
void patch_system::set_ghosted_gridfn_to_constant(fp a, int ghosted_dst_gfn)
{
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	patch& p = ith_patch(pn);
		for (int irho = p.ghosted_min_irho() ;
                     irho <= p.ghosted_max_irho() ; ++irho)
		{
		for (int isigma = p.ghosted_min_isigma() ;
		     isigma <= p.ghosted_max_isigma() ;
		     ++isigma)
		{
		p.ghosted_gridfn(ghosted_dst_gfn, irho,isigma) = a;
		}
		}
	}
}

//******************************************************************************

//
// This function copies one (nominal-grid) gridfn to another.
//
void patch_system::gridfn_copy(int src_gfn, int dst_gfn)
{
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	patch& p = ith_patch(pn);
		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		p.gridfn(dst_gfn, irho,isigma) = p.gridfn(src_gfn, irho,isigma);
		}
		}
	}
}

//******************************************************************************

//
// This function copies one ghosted gridfn to another.
//
void patch_system::ghosted_gridfn_copy(int ghosted_src_gfn,
                                       int ghosted_dst_gfn)
{
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	patch& p = ith_patch(pn);
		for (int irho = p.ghosted_min_irho() ;
                     irho <= p.ghosted_max_irho() ; ++irho)
		{
		for (int isigma = p.ghosted_min_isigma() ;
		     isigma <= p.ghosted_max_isigma() ;
		     ++isigma)
		{
		p.ghosted_gridfn(ghosted_dst_gfn, irho,isigma)
                  = p.ghosted_gridfn(ghosted_src_gfn, irho,isigma);
		}
		}
	}
}

//******************************************************************************

//
// This function copies a ghosted gridfn to a nominal gridfn.
//
void patch_system::ghosted_gridfn_copy_to_nominal(int ghosted_src_gfn,
                                                  int         dst_gfn)
{
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	patch& p = ith_patch(pn);
		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ; isigma <= p.max_isigma() ;
		     ++isigma)
		{
		p.gridfn(dst_gfn, irho,isigma)
                  = p.ghosted_gridfn(ghosted_src_gfn, irho,isigma);
		}
		}
	}
}

//******************************************************************************

//
// This function adds a scalar to a (nominal-grid) gridfn.
//
void patch_system::add_to_gridfn(fp delta, int dst_gfn)
{
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	patch& p = ith_patch(pn);
		for (int irho = p.min_irho() ;
		     irho <= p.max_irho() ;
		     ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		p.gridfn(dst_gfn, irho,isigma) += delta;
		}
		}
	}
}

//******************************************************************************

//
// This function adds a scalar to a ghosted gridfn.
//
void patch_system::add_to_ghosted_gridfn(fp delta, int ghosted_dst_gfn)
{
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	patch& p = ith_patch(pn);
		for (int irho = p.ghosted_min_irho() ;
		     irho <= p.ghosted_max_irho() ;
		     ++irho)
		{
		for (int isigma = p.ghosted_min_isigma() ;
		     isigma <= p.ghosted_max_isigma() ;
		     ++isigma)
		{
		p.ghosted_gridfn(ghosted_dst_gfn, irho,isigma) += delta;
		}
		}
	}
}

//******************************************************************************

//
// This function scales a ghosted gridfn with a scalar.
//
void patch_system::scale_ghosted_gridfn(fp factor, int ghosted_dst_gfn)
{
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	patch& p = ith_patch(pn);
		for (int irho = p.ghosted_min_irho() ;
		     irho <= p.ghosted_max_irho() ;
		     ++irho)
		{
		for (int isigma = p.ghosted_min_isigma() ;
		     isigma <= p.ghosted_max_isigma() ;
		     ++isigma)
		{
		p.ghosted_gridfn(ghosted_dst_gfn, irho,isigma) *= factor;
		}
		}
	}
}

//******************************************************************************

//
// This function computes norms of a nominal-grid gridfn.
//
void patch_system::gridfn_norms(int src_gfn, jtutil::norm<fp>& norms)
	const
{
if (! is_valid_gfn(src_gfn))
   then error_exit(ERROR_EXIT,
"***** patch_system::gridfn_norms(): invalid src_gfn=%d!\n",
		   src_gfn);					/*NOTREACHED*/

norms.reset();

	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	const patch& p = ith_patch(pn);
		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		norms.data(p.gridfn(src_gfn, irho,isigma));
		}
		}
	}
}

//******************************************************************************

//
// This function computes norms of a ghosted-grid gridfn over the
// nominal grid.
//
void patch_system::ghosted_gridfn_norms(int ghosted_src_gfn,
					jtutil::norm<fp>& norms)
	const
{
if (! is_valid_ghosted_gfn(ghosted_src_gfn))
   then error_exit(ERROR_EXIT,
"***** patch_system::gridfn_norms(): invalid ghosted_src_gfn=%d!\n",
		   ghosted_src_gfn);				/*NOTREACHED*/
norms.reset();

	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	const patch& p = ith_patch(pn);
		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		norms.data(p.ghosted_gridfn(ghosted_src_gfn, irho,isigma));
		}
		}
	}
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function computes an approximation to the circumference of a
// surface in the xy, xz, or yz plane.  Note that we compute the full
// circumference all around the sphere, even if the patch system only
// covers a proper subset of this.
//
// We assume that adjacent patches are butt-joined, i.e. that their
// nominal boundaries just touch.
//
// Arguments:
// plane[] = (in) "xy", "xz", or "yz" to specify the integration plane.
// ghosted_radius_gfn = (in) The surface radius.
// g_{xx,xy,xz,yy,yz,zz}_gfn = (in) The xyz 3-metric components.
// method = (in) Selects the integration scheme.
//
fp patch_system::circumference(const char plane[],
			       int ghosted_radius_gfn,
			       int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
					     int g_yy_gfn, int g_yz_gfn,
							   int g_zz_gfn,
			       enum patch::integration_method method)
	const
{
//
// compute arc length around the patch system
//
fp arc_length = 0.0;
    for (int pn = 0 ; pn < N_patches() ; ++pn)
    {
    const patch& p = ith_patch(pn);
    if ((p.ctype() == plane[0]) || (p.ctype() == plane[1]))
       then arc_length += p.plane_arc_length(plane,
					     ghosted_radius_gfn,
					     g_xx_gfn, g_xy_gfn, g_xz_gfn,
						       g_yy_gfn, g_yz_gfn,
								 g_zz_gfn,
					     method);
    }

//
// correct the arc length
// for the fact that the patch system may not cover the full 2-sphere
//
switch	(type())
	{
case patch_system__full_sphere:
	break;
case patch_system__plus_z_hemisphere:
	arc_length *= ((plane[0] == 'x') && (plane[1] == 'y')) ? 1.0 : 2.0;
	break;
case patch_system__plus_xy_quadrant_mirrored:
case patch_system__plus_xy_quadrant_rotating:
	arc_length *= ((plane[0] == 'x') && (plane[1] == 'y')) ? 4.0 : 2.0;
	break;
case patch_system__plus_xz_quadrant_mirrored:
case patch_system__plus_xz_quadrant_rotating:
	arc_length *= ((plane[0] == 'x') && (plane[1] == 'z')) ? 4.0 : 2.0;
	break;
case patch_system__plus_xyz_octant_mirrored:
case patch_system__plus_xyz_octant_rotating:
	arc_length *= 4.0;
	break;
default:
	error_exit(PANIC_EXIT,
"***** patch_system::circumference(): unknown patch system type()=(int)%d!\n"
"                                     (this should never happen!)\n",
		   int(type()));				/*NOTREACHED*/
	}

return arc_length;
}

//******************************************************************************

//
// This function computes an approximation to the (surface) integral of
// a gridfn over the 2-sphere
//	$\int f(\rho,\sigma) \, dA$
//		= \int f(\rho,\sigma) \sqrt{|J|} \, d\rho \, d\sigma$
// where $J$ is the Jacobian of $(x,y,z)$ with respect to $(rho,sigma).
//
// We assume that adjacent patches are butt-joined, i.e. that their
// nominal boundaries just touch.
//
// Arguments:
// unknown_src_gfn = (in) The gridfn to be integrated.  This may be
//			  either nominal-grid or ghosted-grid.
// src_gfn_is_even_across_{xy,xz,yz}_plane
//	= (in) Boolean flags specifying whether the gridfn to be integrated
//	       is even (true) or odd (false) across the corresponding planes.
//	       Only the flags corresponding to boundaries of the patch system
//	       are used.  For example, for a  plus_z_hemisphere  patch system,
//	       only the  src_gfn_is_even_across_xy_plane  flag is used.
// ghosted_radius_gfn = (in) The surface radius.
// g_{xx,xy,xz,yy,yz,zz}_gfn = (in) The xyz 3-metric components.
// method = (in) Selects the integration scheme.
//
fp patch_system::integrate_gridfn(int unknown_src_gfn,
				  bool src_gfn_is_even_across_xy_plane,
				  bool src_gfn_is_even_across_xz_plane,
				  bool src_gfn_is_even_across_yz_plane,
				  int ghosted_radius_gfn,
				  int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
						int g_yy_gfn, int g_yz_gfn,
							      int g_zz_gfn,
				  enum patch::integration_method method)
	const
{
//
// compute integral over patch system
//
fp integral = 0.0;
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	const patch& p = ith_patch(pn);
	integral += p.integrate_gridfn(unknown_src_gfn,
				       ghosted_radius_gfn,
				       g_xx_gfn, g_xy_gfn, g_xz_gfn,
						 g_yy_gfn, g_yz_gfn,
							   g_zz_gfn,
				       method);
	}

//
// correct the integral
// for the fact that the patch system may not cover the full 2-sphere
//
switch	(type())
	{
case patch_system__full_sphere:
	break;
case patch_system__plus_z_hemisphere:
	integral *= src_gfn_is_even_across_xy_plane ? 2.0 : 0.0;
	break;
case patch_system__plus_xy_quadrant_mirrored:
case patch_system__plus_xy_quadrant_rotating:
	integral *= src_gfn_is_even_across_xz_plane ? 2.0 : 0.0;
	integral *= src_gfn_is_even_across_yz_plane ? 2.0 : 0.0;
	break;
case patch_system__plus_xz_quadrant_mirrored:
case patch_system__plus_xz_quadrant_rotating:
	integral *= src_gfn_is_even_across_xy_plane ? 2.0 : 0.0;
	integral *= src_gfn_is_even_across_yz_plane ? 2.0 : 0.0;
	break;
case patch_system__plus_xyz_octant_mirrored:
case patch_system__plus_xyz_octant_rotating:
	integral *= src_gfn_is_even_across_xy_plane ? 2.0 : 0.0;
	integral *= src_gfn_is_even_across_xz_plane ? 2.0 : 0.0;
	integral *= src_gfn_is_even_across_yz_plane ? 2.0 : 0.0;
	break;
default:
	error_exit(PANIC_EXIT,
"***** patch_system::integrate_gridfn(): bad patch system type()=(int)%d!\n"
"                                        (this should never happen!)\n",
		   int(type()));				/*NOTREACHED*/
	}

return integral;
}

fp patch_system::integrate_gridpoint(int unknown_src_gfn,
                                     int ghosted_radius_gfn,
                                     int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
                                                   int g_yy_gfn, int g_yz_gfn,
                                                                 int g_zz_gfn,
                                     enum patch::integration_method method,
                                     int pn, int irho, int isigma)
	const
{
const patch& p = ith_patch(pn);
const fp val = p.integrate_gridpoint(unknown_src_gfn,
                                     ghosted_radius_gfn,
                                     g_xx_gfn, g_xy_gfn, g_xz_gfn,
                                               g_yy_gfn, g_yz_gfn,
                                                         g_zz_gfn,
                                     method,
                                     irho, isigma);

return val;
}
fp patch_system::integrate_correction(bool src_gfn_is_even_across_xy_plane,
                                      bool src_gfn_is_even_across_xz_plane,
                                      bool src_gfn_is_even_across_yz_plane)
	const
{
fp integral = 1.0;
//
// correct the integral
// for the fact that the patch system may not cover the full 2-sphere
//
switch	(type())
	{
case patch_system__full_sphere:
	break;
case patch_system__plus_z_hemisphere:
	integral *= src_gfn_is_even_across_xy_plane ? 2.0 : 0.0;
	break;
case patch_system__plus_xy_quadrant_mirrored:
case patch_system__plus_xy_quadrant_rotating:
	integral *= src_gfn_is_even_across_xz_plane ? 2.0 : 0.0;
	integral *= src_gfn_is_even_across_yz_plane ? 2.0 : 0.0;
	break;
case patch_system__plus_xz_quadrant_mirrored:
case patch_system__plus_xz_quadrant_rotating:
	integral *= src_gfn_is_even_across_xy_plane ? 2.0 : 0.0;
	integral *= src_gfn_is_even_across_yz_plane ? 2.0 : 0.0;
	break;
case patch_system__plus_xyz_octant_mirrored:
case patch_system__plus_xyz_octant_rotating:
	integral *= src_gfn_is_even_across_xy_plane ? 2.0 : 0.0;
	integral *= src_gfn_is_even_across_xz_plane ? 2.0 : 0.0;
	integral *= src_gfn_is_even_across_yz_plane ? 2.0 : 0.0;
	break;
default:
	error_exit(PANIC_EXIT,
"***** patch_system::integrate_gridfn(): bad patch system type()=(int)%d!\n"
"                                        (this should never happen!)\n",
		   int(type()));				/*NOTREACHED*/
	}

return integral;
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function finds what patch contains (the ray from the origin to)
// a given local (x,y,z) position.
//
// If there are multiple patches containing the position, we return the
// one which would still contain it if patches didn't overlap; if multiple
// patches satisfy this criterion then it's arbitrary which one we return.
//
// If no patch contains the position (this can only if the point as at
// the local coordinate origin, or for a non--full-sphere patch system),
// then we return a NULL pointer.
//
// Arguments:
// (x,y,z) = The local coordinates to be converted.
//
// Results:
// This function returns a reference to the containing patch.
//
const patch* patch_system::patch_containing_local_xyz(fp x, fp y, fp z)
	const
{
if ((x == 0.0) && (y == 0.0) && (z == 0.0))
   then return NULL;

// to which axis is (x,y,z) closest?
// ... or equivalently, which of |x|, |y|, and |z| is largest?
const fp abs_x = jtutil::abs(x);
const fp abs_y = jtutil::abs(y);
const fp abs_z = jtutil::abs(z);

if	((abs_z >= abs_x) && (abs_z >= abs_y))
   then return &plus_or_minus_xyz_patch(z > 0.0, 'z');	// +/- z patch
else if ((abs_x >= abs_y) && (abs_x >= abs_z))
   then return &plus_or_minus_xyz_patch(x > 0.0, 'x');	// +/- x patch
else if ((abs_y >= abs_x) && (abs_y >= abs_z))
   then return &plus_or_minus_xyz_patch(y > 0.0, 'y');	// +/- y patch
else	error_exit(ERROR_EXIT,
"***** patch_system::patch_containing_local_xyz():\n"
"        unknown (wierd!) ordering of |x|, |y|, and |z|!\n"
"        (this should never happen!)\n"
"        [local] (x,y,z)=(%g,%g,%g)\n"
,
		   double(x), double(y), double(z));		/*NOTREACHED*/
}

//******************************************************************************

//
// This function computes the radius of a patch-system 2-surface in the
// direction of a specified local (x,y,z) point, taking into account any
// patch-system symmetries.  If the point coincides with the local origin,
// we return the dummy value 1.0.
//
// Due to the surface-interpolator overhead, repeatedly calling this
// function is rather inefficient.
//
fp patch_system::radius_in_local_xyz_direction(int ghosted_radius_gfn,
					       fp x, fp y, fp z)
	const
{
if ((x == 0.0) && (y == 0.0) && (z == 0.0))
   then return 1.0;

//
// apply symmetries to map (x,y,z) into that part of the 2-sphere
// which is covered by the patch system
//
switch	(type())
	{
case patch_system__full_sphere:
	break;
case patch_system__plus_z_hemisphere:
	z = fabs(z);
	break;
case patch_system__plus_xy_quadrant_mirrored:
case patch_system__plus_xy_quadrant_rotating:
	x = fabs(x);
	y = fabs(y);
	break;
case patch_system__plus_xz_quadrant_mirrored:
case patch_system__plus_xz_quadrant_rotating:
	x = fabs(x);
	z = fabs(z);
	break;
case patch_system__plus_xyz_octant_mirrored:
case patch_system__plus_xyz_octant_rotating:
	x = fabs(x);
	y = fabs(y);
	z = fabs(z);
	break;
default:
	error_exit(PANIC_EXIT,
"***** patch_system::radius_in_local_xyz_direction():\n"
"        unknown patch system type()=(int)%d!\n"
"        (this should never happen!)\n",
		   int(type()));				/*NOTREACHED*/
	}

const patch* p_ptr = patch_containing_local_xyz(x, y, z);
if (p_ptr == NULL)
   then error_exit(ERROR_EXIT,
"***** patch_system::radius_in_local_xyz_direction():\n"
"                    can't find containing patch!\n"
"                    (this should never happen!)\n"
"                    [local] (x,y,z)=(%g,%g,%g)\n"
		   ,
		   double(x), double(y), double(z));		/*NOTREACHED*/

const patch& p = *p_ptr;
const fp   rho = p.  rho_of_xyz(x,y,z);
const fp sigma = p.sigma_of_xyz(x,y,z);


//
// Set up the surface interpolator to interpolate the surface radius
// gridfn to the (rho,sigma) coordinates:
//
// Notes on the interpolator setup:
// * The interpolator assumes Fortran subscripting, so we take the
//   coordinates in the order (sigma,rho) to match our C subscripting
//   in the patch system.
// * To avoid having to set up min/max array subscripts in the parameter
//   table, we treat the patch as using 0-origin (integer) array subscripts
//   (irho - ghosted_min_irho(), isigma - ghosted_min_isigma()).  However,
//   we use the usual floating-point coordinates.
//

const int N_dims = 2;
const CCTK_REAL coord_origin[N_dims]
	= { p.ghosted_min_sigma(), p.ghosted_min_rho() };
const CCTK_REAL coord_delta[N_dims]
	= { p.delta_sigma(), p.delta_rho() };

const int N_interp_points = 1;
const int interp_coords_type_code = CCTK_VARIABLE_REAL;
const void* const interp_coords[N_dims]
	= { static_cast<const void*>(&sigma), static_cast<const void*>(&rho) };

const int N_input_arrays = 1;
const CCTK_INT input_array_dims[N_dims]
	= { p.ghosted_N_isigma(), p.ghosted_N_irho() };
const CCTK_INT input_array_type_codes[N_input_arrays] = { CCTK_VARIABLE_REAL };
const void* const input_arrays[N_input_arrays]
	= {
	  static_cast<const void*>(
	    p.ghosted_gridfn_data_array(ghosted_radius_gfn)
				  )
	  };

const int N_output_arrays = 1;
const CCTK_INT output_array_type_codes[N_output_arrays]
	= { CCTK_VARIABLE_REAL };
fp xyz_radius;
void* const output_arrays[N_output_arrays]
	= { static_cast<void*>(&xyz_radius) };

const int status
  = AHFD::CCTK_InterpLocalUniform(N_dims,
				  surface_interp_handle_,
				  surface_interp_par_table_handle_,
				  coord_origin, coord_delta,
				  N_interp_points, interp_coords_type_code,
						   interp_coords,
				  N_input_arrays, input_array_dims,
						  input_array_type_codes,
						  input_arrays,
				  N_output_arrays, output_array_type_codes,
						   output_arrays);
if (status < 0)
   then error_exit(ERROR_EXIT,
"***** patch_system::radius_in_local_xyz_direction():\n"
"                    error return (status=%d) from surface interpolator!\n"
"                    [local] (x,y,z)=(%g,%g,%g)\n"
"                    %s patch (rho,sigma)=(%g,%g)\n"
"                             (drho,dsigma)=(%g,%g)\n"
		   ,
		   status,
		   double(x), double(y), double(z),
		   p.name(), double(rho), double(sigma),
			     double(p.drho_of_rho(rho)),
			     double(p.dsigma_of_sigma(sigma)));	/*NOTREACHED*/

return xyz_radius;
}

struct perpatch {
  std::vector<int> n;
  std::vector<fp> rho, sigma;
  std::vector<fp> radii;
  void reserve (int const npoints)
  {
    n.reserve (npoints);
    rho.reserve (npoints);
    sigma.reserve (npoints);
    radii.reserve (npoints);
  }
  void addpoint (int const n_, fp const rho_, fp const sigma_)
  {
    n    .push_back (n_);
    rho  .push_back (rho_);
    sigma.push_back (sigma_);
    radii.push_back (fp (0.0));
  }
};

void patch_system::radii_in_local_xyz_directions(int const ghosted_radius_gfn,
                                                 int const npoints,
                                                 fp const * const xp,
                                                 fp const * const yp,
                                                 fp const * const zp,
                                                 fp * const radii)
	const
{
  std::vector<perpatch> patchpoints (N_patches());
  for (int pn=0; pn<N_patches(); ++pn) {
    patchpoints.at(pn).reserve (npoints);
  }
  
  for (int n=0; n<npoints; ++n) {
    fp x = xp[n];
    fp y = yp[n];
    fp z = zp[n];
    
    if ((x == 0.0) && (y == 0.0) && (z == 0.0)) {
      radii[n] = 1.0;
    } else {
      //
      // apply symmetries to map (x,y,z) into that part of the
      // 2-sphere which is covered by the patch system
      //
      switch (type())
      {
      case patch_system__full_sphere:
	break;
      case patch_system__plus_z_hemisphere:
	z = fabs(z);
	break;
      case patch_system__plus_xy_quadrant_mirrored:
      case patch_system__plus_xy_quadrant_rotating:
	x = fabs(x);
	y = fabs(y);
	break;
      case patch_system__plus_xz_quadrant_mirrored:
      case patch_system__plus_xz_quadrant_rotating:
	x = fabs(x);
	z = fabs(z);
	break;
      case patch_system__plus_xyz_octant_mirrored:
      case patch_system__plus_xyz_octant_rotating:
	x = fabs(x);
	y = fabs(y);
	z = fabs(z);
	break;
      default:
	error_exit(PANIC_EXIT,
"***** patch_system::radii_in_local_xyz_directions():\n"
"        unknown patch system type()=(int)%d!\n"
"        (this should never happen!)\n",
		   int(type()));				/*NOTREACHED*/
      }
      
      const patch* p_ptr = patch_containing_local_xyz(x, y, z);
      if (p_ptr == NULL)
        then error_exit(ERROR_EXIT,
"***** patch_system::radius_in_local_xyz_direction():\n"
"                    can't find containing patch!\n"
"                    (this should never happen!)\n"
"                    [local] (x,y,z)=(%g,%g,%g)\n"
		   ,
		   double(x), double(y), double(z));		/*NOTREACHED*/
      
      const patch& p = *p_ptr;
      const fp   rho = p.  rho_of_xyz(x,y,z);
      const fp sigma = p.sigma_of_xyz(x,y,z);
      
      patchpoints.at(p.patch_number()).addpoint (n, rho, sigma);
    }
  } // for n
  
  for (int pn=0; pn<N_patches(); ++pn) {
    const patch& p = ith_patch(pn);
    perpatch& localpatchpoints = patchpoints.at(pn);
    
    // Set up the surface interpolator to interpolate the surface
    // radius gridfn to the (rho,sigma) coordinates:
    //
    // Notes on the interpolator setup:
    // * The interpolator assumes Fortran subscripting, so we take the
    //   coordinates in the order (sigma,rho) to match our C
    //   subscripting in the patch system.
    // * To avoid having to set up min/max array subscripts in the
    //   parameter table, we treat the patch as using 0-origin
    //   (integer) array subscripts (irho - ghosted_min_irho(), isigma
    //   - ghosted_min_isigma()).  However, we use the usual
    //   floating-point coordinates.
    //
    
    const int N_dims = 2;
    const CCTK_REAL coord_origin[N_dims]
      = { p.ghosted_min_sigma(), p.ghosted_min_rho() };
    const CCTK_REAL coord_delta[N_dims]
      = { p.delta_sigma(), p.delta_rho() };
    
    const int N_interp_points = localpatchpoints.radii.size();
    const int interp_coords_type_code = CCTK_VARIABLE_REAL;
    const void* const interp_coords[N_dims]
      = { static_cast<const void*>(& localpatchpoints.sigma.front()),
          static_cast<const void*>(& localpatchpoints.rho  .front()) };
    
    const int N_input_arrays = 1;
    const CCTK_INT input_array_dims[N_dims]
      = { p.ghosted_N_isigma(), p.ghosted_N_irho() };
    const CCTK_INT input_array_type_codes[N_input_arrays]
      = { CCTK_VARIABLE_REAL };
    const void* const input_arrays[N_input_arrays]
      = { static_cast<const void*>
          (p.ghosted_gridfn_data_array(ghosted_radius_gfn)) };
    
    const int N_output_arrays = 1;
    const CCTK_INT output_array_type_codes[N_output_arrays]
      = { CCTK_VARIABLE_REAL };
    void* const output_arrays[N_output_arrays]
      = { static_cast<void*>(& localpatchpoints.radii.front()) };
    
    const int status = AHFD::CCTK_InterpLocalUniform
      (N_dims,
       surface_interp_handle_,
       surface_interp_par_table_handle_,
       coord_origin, coord_delta,
       N_interp_points,
       interp_coords_type_code, interp_coords,
       N_input_arrays,
       input_array_dims, input_array_type_codes, input_arrays,
       N_output_arrays,
       output_array_type_codes, output_arrays);
    if (status < 0)
      then error_exit(ERROR_EXIT,
"***** patch_system::radii_in_local_xyz_directions():\n"
"                    error return (status=%d) from surface interpolator!\n"
		   ,
                      status);	/*NOTREACHED*/
    
    for (size_t nn=0; nn<localpatchpoints.radii.size(); ++nn) {
      radii[localpatchpoints.n.at(nn)] = localpatchpoints.radii.at(nn);
    }
    
  } // for p
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function prints an unknown-grid gridfn in ASCII format to a
// named output file.  The output format is suitable for a gnuplot
// 'splot' command.  (Individual patches may be selected with the
//  select.patch  program (perl script).)  The output format is either
//	# print_xyz_flag == false
//	dpx	dpy	gridfn
// or
//	# print_xyz_flag == true
//	dpx	dpy	gridfn   global_x   global_y   global_z
// where global_[xyz] are derived from the angular position and a
// specified (unknown-grid) radius gridfn.
//
void patch_system::print_unknown_gridfn
	(bool ghosted_flag, int unknown_gfn,
	 bool print_xyz_flag, bool radius_is_ghosted_flag,
			      int unknown_radius_gfn,
	 const char output_file_name[], bool want_ghost_zones)
	const
{
if (want_ghost_zones && !ghosted_flag)
   then error_exit(PANIC_EXIT,
"***** patch_system::print_unknown_gridfn(unknown_gfn=%d):\n"
"        can't have want_ghost_zones && !ghosted_flag !\n"
,
		   unknown_gfn);				/*NOTREACHED*/
if (want_ghost_zones && print_xyz_flag && !radius_is_ghosted_flag)
   then error_exit(PANIC_EXIT,
"***** patch_system::print_unknown_gridfn(unknown_gfn=%d):\n"
"        can't have want_ghost_zones && print_xyz_flag\n"
"        && !radius_is_ghosted_flag!\n"
"        unknown_radius_gfn=%d\n"
,
		   unknown_gfn,
		   unknown_radius_gfn);				/*NOTREACHED*/

FILE *output_fp = fopen(output_file_name, "w");
if (output_fp == NULL)
   then error_exit(ERROR_EXIT,
"***** patch_system::print_unknown_gridfn(unknown_gfn=%d):\n"
"        can't open output file \"%s\"\n!"
,
		   unknown_gfn,
		   output_file_name);				/*NOTREACHED*/

fprintf(output_fp, "# N_patches = %d\n", N_patches());
fprintf(output_fp, "# origin = %.15g %.15g %.15g\n",
        double(origin_x()), double(origin_y()), double(origin_z()));
fprintf(output_fp, "\n");

	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	const patch& p = ith_patch(pn);

	fprintf(output_fp, "### %s patch\n", p.name());
	fprintf(output_fp, "# N_rho = %d\n",
		p.effective_N_irho(want_ghost_zones));
	fprintf(output_fp, "# N_sigma = %d\n",
		p.effective_N_isigma(want_ghost_zones));
	fprintf(output_fp, "# %s_gfn=%d\n",
		(ghosted_flag ? "ghosted" : "nominal"), unknown_gfn);
	fprintf(output_fp, "# dpx = %s\n", p.name_of_dpx());
	fprintf(output_fp, "# dpy = %s\n", p.name_of_dpy());
	fprintf(output_fp, "#\n");
	fprintf(output_fp,
		print_xyz_flag
		? "# dpx\tdpy\tgridfn\tglobal_x\tglobal_y\tglobal_z\n"
		: "# dpx\tdpy\tgridfn\n");

		for (int irho = p.effective_min_irho(want_ghost_zones) ;
		     irho <= p.effective_max_irho(want_ghost_zones) ;
		     ++irho)
		{
		for (int isigma = p.effective_min_isigma(want_ghost_zones) ;
		     isigma <= p.effective_max_isigma(want_ghost_zones) ;
		     ++isigma)
		{
		const fp rho = p.rho_of_irho(irho);
		const fp sigma = p.sigma_of_isigma(isigma);
		const fp dpx = p.dpx_of_rho_sigma(rho, sigma);
		const fp dpy = p.dpy_of_rho_sigma(rho, sigma);
		fprintf(output_fp,
			"%g\t%g\t%#.15g",
			double(dpx), double(dpy),
                        double(p.unknown_gridfn(ghosted_flag,
                                                unknown_gfn, irho,isigma)));
		if (print_xyz_flag)
		   then {
			const fp r = p.unknown_gridfn(radius_is_ghosted_flag,
						      unknown_radius_gfn,
						      irho,isigma);
			fp local_x, local_y, local_z;
			p.xyz_of_r_rho_sigma(r,rho,sigma,
					     local_x,local_y,local_z);
			const fp global_x = origin_x() + local_x;
			const fp global_y = origin_y() + local_y;
			const fp global_z = origin_z() + local_z;
			fprintf(output_fp,
				"\t%#.10g\t%#.10g\t%#.10g",
				double(global_x), double(global_y),
                                double(global_z));
			}
		fprintf(output_fp, "\n");
		}
		fprintf(output_fp, "\n");
		}
	fprintf(output_fp, "\n");
	}

fclose(output_fp);
}

//******************************************************************************

//
// This function reads an unknown-grid gridfn in ASCII format from
// a named input file.  Comments ('#' in column 1) and blank lines
// are ignored, otherwise the input format matches that written by
// print_unknown_gridfn(): the first 3 numbers on each line are taken
// to be dpx, dpy, and the gridfn value; anything else on the line is
// ignored.
//
void patch_system::read_unknown_gridfn(bool ghosted_flag, int unknown_gfn,
				       const char input_file_name[],
				       bool want_ghost_zones)
{
if (want_ghost_zones && !ghosted_flag)
   then error_exit(PANIC_EXIT,
"***** patch_system::read_unknown_gridfn(unknown_gfn=%d):\n"
"        can't have want_ghost_zones && !ghosted_flag !\n"
,
		   unknown_gfn);				/*NOTREACHED*/

FILE *input_fp = fopen(input_file_name, "r");
if (input_fp == NULL)
   then error_exit(ERROR_EXIT,
"***** patch_system::read_unknown_gridfn(unknown_gfn=%d):\n"
"        can't open input file \"%s\"\n!"
,
		   unknown_gfn,
		   input_file_name);				/*NOTREACHED*/

int line_number = 1;
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	patch& p = ith_patch(pn);

	for (int irho = p.effective_min_irho(want_ghost_zones) ;
	     irho <= p.effective_max_irho(want_ghost_zones) ;
	     ++irho)
	{
	for (int isigma = p.effective_min_isigma(want_ghost_zones) ;
	     isigma <= p.effective_max_isigma(want_ghost_zones) ;
	     ++isigma)
	{
	const fp rho   = p.rho_of_irho(irho);
	const fp sigma = p.sigma_of_isigma(isigma);
	const fp dpx = p.dpx_of_rho_sigma(rho, sigma);
	const fp dpy = p.dpy_of_rho_sigma(rho, sigma);

	const int buffer_size = 250;
	char buffer[buffer_size];
	// read/discard comments and blank lines
		do
		{
		if (fgets(buffer, buffer_size, input_fp) == NULL)
		   then error_exit(ERROR_EXIT,
"***** patch::read_unknown_gridfn(%s patch, unknown_gfn=%d):\n"
"        I/O error or unexpected end-of-file on input!\n"
"        at irho=%d of [%d,%d], isigma=%d of [%d,%d]\n"
"        dpx=%g dpy=%g\n"
,
				   p.name(), unknown_gfn,
				   int(irho),
                                      p.effective_min_irho(want_ghost_zones),
				      p.effective_max_irho(want_ghost_zones),
				   int(isigma),
				      p.effective_min_isigma(want_ghost_zones),
				      p.effective_max_isigma(want_ghost_zones),
				   double(dpx), double(dpy));	/*NOTREACHED*/
		++line_number;
		} while ((buffer[0] == '#') || (buffer[0] == '\n'));

	double read_dpx, read_dpy, read_gridfn_value;
	if (sscanf(buffer, "%lf %lf %lf",
		   &read_dpx, &read_dpy, &read_gridfn_value) != 3)
	   then error_exit(ERROR_EXIT,
"***** patch::read_unknown_gridfn(%s patch, unknown_gfn=%d):\n"
"        bad input data at input line %d!\n"
,
			   p.name(), unknown_gfn,
			   line_number);			/*NOTREACHED*/
	if (! (    jtutil::fuzzy<fp>::EQ(read_dpx,dpx)
		&& jtutil::fuzzy<fp>::EQ(read_dpy,dpy)    ) )
	   then error_exit(ERROR_EXIT,
"***** patch::read_unknown_gridfn(%s patch, unknown_gfn=%d):\n"
"        wrong (dpx,dpy) at input line %d!\n"
"        expected (%g,%g)\n"
"        read     (%g,%g)\n"
,
			   p.name(), unknown_gfn,
			   line_number,
			   double(dpx), double(dpy),
			   double(read_dpx), double(read_dpy));	/*NOTREACHED*/

	p.unknown_gridfn(ghosted_flag,
			 unknown_gfn, irho,isigma) = read_gridfn_value;
	}
	}

	}

fclose(input_fp);
}

//******************************************************************************

#ifdef HAVE_CAPABILITY_HDF5

static void WriteAttribute (const hid_t dataset, const char* name, char value);
static void WriteAttribute (const hid_t dataset, const char* name, const char* values);
static void WriteAttribute (const hid_t dataset, const char* name, const char* values, int nvalues);

static void WriteAttribute (const hid_t dataset, const char* const name, const char value)
{
  WriteAttribute (dataset, name, &value, 1);
}

static void WriteAttribute (const hid_t dataset, const char* const name, const char* const values)
{
  WriteAttribute (dataset, name, values, strlen(values));
}

static void WriteAttribute (const hid_t dataset, const char* const name, const char* const values, const int nvalues)
{
  assert (dataset>=0);
  assert (name);
  assert (values);
  assert (nvalues>=0);
  
  herr_t herr;
  
  const hid_t dataspace = H5Screate (H5S_SCALAR);
  assert (dataspace>=0);
  
  const hid_t datatype = H5Tcopy (H5T_C_S1);
  assert (datatype>=0);
  herr = H5Tset_size (datatype, nvalues);
  assert (!herr);
  
  const hid_t attribute = H5Acreate (dataset, name, datatype, dataspace, H5P_DEFAULT);
  assert (attribute>=0);
  herr = H5Awrite (attribute, datatype, values);
  assert (!herr);
  herr = H5Aclose (attribute);
  assert (!herr);
  
  herr = H5Tclose (datatype);
  assert (!herr);
  
  herr = H5Sclose (dataspace);
  assert (!herr);
}

#endif

//
// This function output an unknown-grid gridfn in HDF5 format to a
// named output file.
//
void patch_system::output_unknown_gridfn
	(bool ghosted_flag, int unknown_gfn,
         const char gfn_name[],
	 bool output_xyz_flag, bool radius_is_ghosted_flag,
                               int unknown_radius_gfn,
	 const char output_file_name[], bool want_ghost_zones)
	const
{
if (want_ghost_zones && !ghosted_flag)
   then error_exit(PANIC_EXIT,
"***** patch_system::output_unknown_gridfn(unknown_gfn=%d):\n"
"        can't have want_ghost_zones && !ghosted_flag !\n"
,
		   unknown_gfn);				/*NOTREACHED*/
if (want_ghost_zones && output_xyz_flag && !radius_is_ghosted_flag)
   then error_exit(PANIC_EXIT,
"***** patch_system::output_unknown_gridfn(unknown_gfn=%d):\n"
"        can't have want_ghost_zones && output_xyz_flag\n"
"        && !radius_is_ghosted_flag!\n"
"        unknown_radius_gfn=%d\n"
,
		   unknown_gfn,
		   unknown_radius_gfn);				/*NOTREACHED*/

#ifdef HAVE_CAPABILITY_HDF5
  
  using std::ostringstream;
  using std::string;
  using std::vector;
    
  herr_t herr;
  
  hid_t writer;

  // TODOMARKS
  // need to be careful here
  // Get grid hierarchy extentsion from IOUtil
  //const ioGH * const iogh = (const ioGH *)CCTK_GHExtension (cctkGH, "IO");
  //assert (iogh);
  
  struct stat fileinfo;
  writer = H5Fcreate (output_file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  // if (! iogh->recovered || stat(output_file_name, &fileinfo) != 0) {
  //   writer = H5Fcreate (output_file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  // } else {
  //   writer = H5Fopen (output_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
  // }
  assert (writer >= 0);
  
  hid_t datatype;
  if (sizeof(fp) == sizeof(float)) {
    datatype = H5T_NATIVE_FLOAT;
  } else if (sizeof(fp) == sizeof(double)) {
    datatype = H5T_NATIVE_DOUBLE;
  } else if (sizeof(fp) == sizeof(long double)) {
    datatype = H5T_NATIVE_LDOUBLE;
  } else {
    assert (0);
  }
  assert (datatype >= 0);
  
  for (int pn=0; pn<N_patches(); ++pn) {
    const patch& p = ith_patch(pn);
    
    hsize_t shape[2];
    shape[0] = p.effective_N_isigma(want_ghost_zones);
    shape[1] = p.effective_N_irho(want_ghost_zones);
    
    hid_t dataspace = H5Screate_simple (2, shape, NULL);
    assert (dataspace >= 0);
    
    
    // TODOMARKS
    // could add step info in the future
    ostringstream datasetnamebuf;
    datasetnamebuf << gfn_name
                   << " it=" << 666
                   << " patch=" << pn;
    string datasetnamestr = datasetnamebuf.str();
    assert (datasetnamestr.size() <= 256); // limit dataset name length
    const char * const datasetname = datasetnamestr.c_str();
    assert (datasetname);
    
    ostringstream datasetnamebufx;
    datasetnamebufx << "x"
                    << " it=" << 666
                    << " patch=" << pn;
    string datasetnamestrx = datasetnamebufx.str();
    assert (datasetnamestrx.size() <= 256); // limit dataset name length
    const char * const datasetnamex = datasetnamestrx.c_str();
    assert (datasetnamex);
    
    ostringstream datasetnamebufy;
    datasetnamebufy << "y"
                    << " it=" << 666
                    << " patch=" << pn;
    string datasetnamestry = datasetnamebufy.str();
    assert (datasetnamestry.size() <= 256); // limit dataset name length
    const char * const datasetnamey = datasetnamestry.c_str();
    assert (datasetnamey);
    
    ostringstream datasetnamebufz;
    datasetnamebufz << "z"
                    << " it=" << 666
                    << " patch=" << pn;
    string datasetnamestrz = datasetnamebufz.str();
    assert (datasetnamestrz.size() <= 256); // limit dataset name length
    const char * const datasetnamez = datasetnamestrz.c_str();
    assert (datasetnamez);
    
    
    
    hid_t dataset = H5Dcreate (writer, datasetname, datatype, dataspace, H5P_DEFAULT);
    assert (dataset >= 0);
    
    hid_t datasetx;
    hid_t datasety;
    hid_t datasetz;
    if (output_xyz_flag) {
      
      datasetx = H5Dcreate (writer, datasetnamex, datatype, dataspace, H5P_DEFAULT);
      assert (datasetx >= 0);
      
      datasety = H5Dcreate (writer, datasetnamey, datatype, dataspace, H5P_DEFAULT);
      assert (datasety >= 0);
      
      datasetz = H5Dcreate (writer, datasetnamez, datatype, dataspace, H5P_DEFAULT);
      assert (datasetz >= 0);
      
    } // if output_xyf_flag
    
    
    
    vector<fp> tmpdata;
    vector<fp> tmpx, tmpy, tmpz;
    tmpdata.resize (shape[0] * shape[1]);
    if (output_xyz_flag) {
      tmpx.resize (shape[0] * shape[1]);
      tmpy.resize (shape[0] * shape[1]);
      tmpz.resize (shape[0] * shape[1]);
    }
    
    int elt = 0;
    for (int irho = p.effective_min_irho(want_ghost_zones);
         irho <= p.effective_max_irho(want_ghost_zones);
         ++irho)
    {
      for (int isigma = p.effective_min_isigma(want_ghost_zones);
           isigma <= p.effective_max_isigma(want_ghost_zones);
           ++isigma)
      {
        
        tmpdata.at(elt) = p.unknown_gridfn (ghosted_flag, unknown_gfn, irho,isigma);
        if (output_xyz_flag) {
          const fp r = p.unknown_gridfn(radius_is_ghosted_flag,
                                        unknown_radius_gfn,
                                        irho,isigma);
          const fp rho   = p.rho_of_irho(irho);
          const fp sigma = p.sigma_of_isigma(isigma);
          fp local_x, local_y, local_z;
          p.xyz_of_r_rho_sigma(r,rho,sigma,
                               local_x,local_y,local_z);
          const fp global_x = origin_x() + local_x;
          const fp global_y = origin_y() + local_y;
          const fp global_z = origin_z() + local_z;
          tmpx.at(elt) = global_x;
          tmpy.at(elt) = global_y;
          tmpz.at(elt) = global_z;
        }
        ++elt;
        
      }
    }
    
    
    
    herr = H5Dwrite (dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmpdata.front());
    assert (! herr);
    
    WriteAttribute (dataset, "name", gfn_name);
    
    herr = H5Dclose (dataset);
    assert (! herr);
    
    
    
    if (output_xyz_flag) {
      
      herr = H5Dwrite (datasetx, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmpx.front());
      assert (! herr);
      
      WriteAttribute (datasetx, "name", "x");
      
      herr = H5Dclose (datasetx);
      assert (! herr);
      
      
      
      herr = H5Dwrite (datasety, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmpy.front());
      assert (! herr);
      
      WriteAttribute (datasety, "name", "y");
      
      herr = H5Dclose (datasety);
      assert (! herr);
      
      
      
      herr = H5Dwrite (datasetz, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmpz.front());
      assert (! herr);
      
      WriteAttribute (datasetz, "name", "z");
      
      herr = H5Dclose (datasetz);
      assert (! herr);
      
    } // if output_xyz_flag
    
    
    
    herr = H5Sclose (dataspace);
    assert (! herr);
    
  } // for pn
  
  
  
  herr = H5Fclose (writer);
  assert (! herr);
  
#else
  error_exit(ERROR_EXIT,
"***** patch_system::output_unknown_gridfn(unknown_gfn=%d):\n"
"        no HDF5 support compiled in.  Cannot write output file \"%s\"\n!"
,
             unknown_gfn,
             output_file_name);					/*NOTREACHED*/

#endif
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function "synchronizes" all ghost zones of all patches, i.e. it
// update the ghost-zone values of the specified gridfns via the appropriate
// sequence of symmetry operations and interpatch interpolations.  This
// process is described in detail in the header comments in "ghost_zone.hh".
//
void patch_system::synchronize(int ghosted_min_gfn_to_sync,
			       int ghosted_max_gfn_to_sync)
{
//
// Phase 1:
// Fill in gridfn data at all the non-corner points of symmetry ghost
// zones, using the symmetries to get this data from its "home patch"
// nominal grids.
//
	  {
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	patch& p = ith_patch(pn);
		// n.b. these loops must use _int_ variables for the loop
		//      to terminate!
		for (int want_min = false ; want_min <= true ; ++want_min)
		{
		for (int want_rho = false ; want_rho <= true ; ++want_rho)
		{
		ghost_zone& gz = p.minmax_ang_ghost_zone(want_min, want_rho);
		if (gz.is_symmetry())
		   then gz.synchronize(ghosted_min_gfn_to_sync,
				       ghosted_max_gfn_to_sync,
				       false,	// want corners?
				       true);	// want non-corner?
		}
		}
	}
	  }

//
// Phase 2:
// Fill in gridfn data in all the interpatch ghost zones, using interpatch
// interpolation from neighboring patches as described above.
//
	  {
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	patch& p = ith_patch(pn);
		// n.b. these loops must use _int_ variables for the loop
		//      to terminate!
		for (int want_min = false ; want_min <= true ; ++want_min)
		{
		for (int want_rho = false ; want_rho <= true ; ++want_rho)
		{
		ghost_zone& gz = p.minmax_ang_ghost_zone(want_min, want_rho);
		if (gz.is_interpatch())
		   then gz.synchronize(ghosted_min_gfn_to_sync,
				       ghosted_max_gfn_to_sync);
		}
		}
	}
	  }

//
// Phase 3:
// Fill in gridfn data at all the corner points of symmetry ghost zones,
// using the symmetries to get this data from its "home patch" nominal
// grids or ghost zones.
//
	  {
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	patch& p = ith_patch(pn);
		// n.b. these loops must use _int_ variables for the loop
		//      to terminate!
		for (int want_min = false ; want_min <= true ; ++want_min)
		{
		for (int want_rho = false ; want_rho <= true ; ++want_rho)
		{
		ghost_zone& gz = p.minmax_ang_ghost_zone(want_min, want_rho);
		if (gz.is_symmetry())
		   then gz.synchronize(ghosted_min_gfn_to_sync,
				       ghosted_max_gfn_to_sync,
				       true,	// want corners?
				       false);	// want non-corner?
		}
		}
	}
	  }
}

//******************************************************************************

//
// This function does any precomputation necessary to compute the Jacobian
// of  synchronize() , taking into account synchronize()'s full 3-phase
// algorithm.  In practice, this means it computes the individual Jacobian
// of each ghost zone, and sets  global_{min,max}_ym_ .
//
void patch_system::compute_synchronize_Jacobian(int ghosted_min_gfn_to_sync,
						int ghosted_max_gfn_to_sync)
	const
{
global_min_ym_ = +INT_MAX;
global_max_ym_ = -INT_MAX;
	for (int pn = 0 ; pn < N_patches() ; ++pn)
	{
	const patch& p = ith_patch(pn);
		// n.b. these loops must use _int_ variables for the loop
		//      to terminate!
		for (int want_min = false ; want_min <= true ; ++want_min)
		{
		for (int want_rho = false ; want_rho <= true ; ++want_rho)
		{
		ghost_zone& gz = p.minmax_ang_ghost_zone(want_min, want_rho);
		gz.compute_Jacobian(ghosted_min_gfn_to_sync,
				    ghosted_max_gfn_to_sync);

		global_min_ym_ = jtutil::min(global_min_ym_,
					     gz.Jacobian_min_y_ipar_m());
		global_max_ym_ = jtutil::max(global_max_ym_,
					     gz.Jacobian_max_y_ipar_m());
		}
		}
	}
}

//******************************************************************************

//
// Given that  compute_synchronize_Jacobian()  has been called, this
// function computes the global min/max  m  over all ghost zone points.
// This is useful for sizing the buffer for synchronize_Jacobian().
//
void patch_system::synchronize_Jacobian_global_minmax_ym
	(int& min_ym, int& max_ym)
	const
{
min_ym = global_min_ym_;
max_ym = global_max_ym_;
}

//******************************************************************************

//
// Given that  compute_synchronize_Jacobian()  has been called, this
// function computes a single row of the Jacobian, taking into account
// synchronize()'s 3-phase algorithm:
// - It returns the edge to which the y point belongs (the caller can get
//   the patch from this edge).
// - It stores y_iperp and y_posn and min/max ym in the named arguments.
// - It stores the Jacobian elements
//	partial synchronize() gridfn(ghosted_gfn, px, x_iperp, x_ipar)
//	-------------------------------------------------------------
//	     partial gridfn(ghosted_gfn, py, y_iperp, y_posn+ym)
//   in the caller-supplied buffer
//	Jacobian_buffer(ym)
//   for each  ym  in the min/max ym range.
//
// In practice, the main task of this function is taking into account
// synchronize()'s 3-phase algorithm.  There are several cases:
// - ghost zone is symmetry && x point is in non-corner
//   ==> x value was computed by a phase 1 symmetry operation,
//       using (only) nominal-grid data
//   ==> overall Jacobian(x,y) = ghost zone Jacobian(x,y)
// - ghost zone is symmetry && x point is in corner
//   --> x value was computed by a phase 3 symmetry operation,
//       from some point (call it z),
//       ==> overall Jacobian(x,y) = overall Jacobian(z,y)
//       ==> call this function recursively to get z's Jacobian
//           (z must be in the noncorner part of some ghost zone,
//	      so this won't lead to infinite recursion)
// - ghost zone is interpatch
//   ==> x value was computed by a phase 2 interpatch interpolation
//       - using (only) nominal-grid data
//         ==> overall Jacobian(x,y) = ghost zone Jacobian(x,y)
//       - using a mixture of nominal-grid data
//         and data computed by a phase 1 symmetry operation
//         ==> overall Jacobian(x,y) = "fold" ghost zone Jacobian(x,y)
//                                     to take the phase 1 symmetry
//                                     operation into account
//
const patch_edge&
  patch_system::synchronize_Jacobian(const ghost_zone& xgz,
				     int x_iperp, int x_ipar,
				     int& y_iperp,
				     int& y_posn, int& min_ym, int& max_ym,
				     jtutil::array1d<fp>& Jacobian_buffer)
	const
{
const patch_edge& xe = xgz.my_edge();

if (xgz.is_symmetry() && xe.ipar_is_in_noncorner(x_ipar))
   then {
	// ghost zone is symmetry && x point is in non-corner
	// ==> x value was computed by a phase 1 symmetry operation,
	//     using (only) nominal-grid data
	// ==> overall Jacobian(x,y) = ghost zone Jacobian(x,y)
	return ghost_zone_Jacobian(xgz,
				   x_iperp, x_ipar,
				   y_iperp,
				   y_posn, min_ym, max_ym,
				   Jacobian_buffer);
	}

else if (xgz.is_symmetry() && xe.ipar_is_in_corner(x_ipar))
   then {
	// ghost zone is symmetry && x point is in corner
	// --> x value was computed by a phase 3 symmetry operation,
	//     from some point (call it z),
	//     ==> overall Jacobian(x,y) = overall Jacobian(z,y)
	//     ==> call this function recursively to get z's Jacobian
	//         (z must be in the noncorner part of some ghost zone,
	//          so this won't lead to infinite recursion)

	const patch&      zp = xgz.other_patch();
	const patch_edge& ze = xgz.other_edge();
	const symmetry_ghost_zone& xsgz = xgz.cast_to_symmetry_ghost_zone();
	const int z_iperp = xsgz.iperp_map_of_iperp(x_iperp);
	const int z_ipar  = xsgz.ipar_map_of_ipar  (x_ipar);

	//
	// Computing z's edge/ghost zone is tricky.  For example:
	//	                        |
	//	           p1         e3|e4          p2
	//	                        |
	//	                        |   z
	//	-----------e1-----------+------------e2---------
	//	                        |   x
	//	                        |
	// Here the point x in the corner of p1's e1 ghost zone,
	// is computed by the phase 3 symmetry operation (a reflection
	// about e1) from z.  Thus zp == p1 and ze == e1.
	//
	// But we need to "turn the corner" to compute z's "true" edge
	// e3 (so we can recursively call this function to compute z's
	// Jacobian).  Thus we explicitly check which ghost zone of p1
	// (here the e3 ghost zone) contains the point z.
	//
	const int z_irho   = ze.  irho_of_iperp_ipar(z_iperp, z_ipar);
	const int z_isigma = ze.isigma_of_iperp_ipar(z_iperp, z_ipar);
	const ghost_zone& true_zgz
		= zp.ghost_zone_containing_noncorner_point(z_irho, z_isigma);
	const patch_edge& true_ze = true_zgz.my_edge();
	const int true_z_iperp = true_ze.iperp_of_irho_isigma(z_irho, z_isigma);
	const int true_z_ipar  = true_ze. ipar_of_irho_isigma(z_irho, z_isigma);

	// make sure we have the right ghost zone!
	assert( true_zgz.is_in_ghost_zone(true_z_iperp, true_z_ipar) );

	return synchronize_Jacobian(true_zgz,
				    true_z_iperp, true_z_ipar,
				    y_iperp,
				    y_posn, min_ym, max_ym,
				    Jacobian_buffer);
	}

else if (xgz.is_interpatch())
   then {
	// ghost zone is interpatch
	// ==> x value was computed by a phase 2 interpatch interpolation
	//     - using (only) nominal-grid data
	//       ==> overall Jacobian(x,y) = ghost zone Jacobian(x,y)
	//     - using a mixture of nominal-grid data
	//       and data computed by a phase 1 symmetry operation
	//       ==> overall Jacobian(x,y) = "fold" ghost zone Jacobian(x,y)
	//                                   to take the phase 1 symmetry
	//                                   operation into account
	//
	// For example,
	//	                        |
	//	           xp         xe|ye a        yp
	//	                        |   b
	//	                        |  xc
	//	----------xse-----------+---d-------yse----------
	//	                        |   e
	//	                        |
	// here point x is computed by interpatch-interpolating in the
	// par direction from the 5 y points abcde.  e is outside the
	// nominal grid, so its Jacobian must be "folded" over to c.
	// Notice that this "folding" must be done about the edge yse,
	// *not* about ye itself.

	// Jacobian of the phase 2 interpatch interpolation
	const patch_edge& ye = ghost_zone_Jacobian(xgz,
						   x_iperp, x_ipar,
						   y_iperp,
						   y_posn, min_ym, max_ym,
						   Jacobian_buffer);
	const int min_y_ipar = y_posn + min_ym;
	const int max_y_ipar = y_posn + max_ym;

	// fold any points in the Jacobian outside the nominal grid
	if (ye.ipar_is_in_min_ipar_corner(min_y_ipar))
	   then {
		fold_Jacobian(ye, ye.min_par_adjacent_edge(),
			      y_iperp,
			      y_posn, min_ym, max_ym,
			      min_ym, ye.min_ipar_corner__max_ipar() - y_posn,
			      Jacobian_buffer);
		min_ym = ye.min_ipar_without_corners() - y_posn;
		}
	if (ye.ipar_is_in_max_ipar_corner(max_y_ipar))
	   then {
		fold_Jacobian(ye, ye.max_par_adjacent_edge(),
			      y_iperp,
			      y_posn, min_ym, max_ym,
			      ye.max_ipar_corner__min_ipar() - y_posn, max_ym,
			      Jacobian_buffer);
		max_ym = ye.max_ipar_without_corners() - y_posn;
		}

	return ye;
	}

else	error_exit(PANIC_EXIT,
"***** patch_system::synchronize_Jacobian():\n"
"        don't know what to do with ghost zone (this should never happen)!\n"
"        xgz.my_patch()=\"%s\"    xe=xgz.my_edge()=\"%s\"\n"
"        xgz.other_patch()=\"%s\"    xgz.other_edge()=\"%s\"\n"
"        xgz.is_symmetry()=(int)%d xgz.is_interpatch()=(int)%d\n"
"        x_iperp=%d x_ipar=%d\n"
"        xe.ipar_is_in_{min,max}_ipar_corner(x_ipar)=(int){%d,%d}\n"
"        xe.ipar_is_in_{corner,noncorner}(x_ipar)=(int){%d,%d}\n"
,
		   xgz.my_patch().name(),    xe.name(),
		   xgz.other_patch().name(), xgz.other_edge().name(),
		   int(xgz.is_symmetry()), int(xgz.is_interpatch()),
		   x_iperp, x_ipar,
		   xe.ipar_is_in_min_ipar_corner(x_ipar),
		   xe.ipar_is_in_max_ipar_corner(x_ipar),
		   xe.ipar_is_in_corner(x_ipar),
		   xe.ipar_is_in_noncorner(x_ipar));		/*NOTREACHED*/
}

//******************************************************************************

//
// This function "folds" part of a(n interpatch) Jacobian row to take
// a symmetry operation into account.  For example:
//	         |
//	         |e_Jac
//	         |             p
//	         |   a
//	         |   b
//	         |   c=y
//	---------+---d-------e_fold-------
//	         |   e=x    sgz_fold
//	         |
// Here the Jacobian abcde is to be "folded", because e is outside the
// nominal grid (its Jacobian must be "folded" over to c).
//
// Notice that the folding (about the edge e_fold) is in the par direction
// with respect to e_Jac, but the perp direction with respect to e_fold.
// Since e_fold and e_Jac are adjacent edges, 
//	e_Jac (iperp,ipar) == e_fold (ipar,iperp)
//
// Arguments:
// e_Jac = edge which the Jacobian lies along
// e_fold = edge about which to fold; the corresponding ghost zone must be
//	    symmetry ghost zone, and at present we only support the case
//	    where this is a "local" (mirror-image) symmetry ghost zone
// iperp = iperp-wrt-e_Jac coordinate of Jacobian
// posn = ipar-wrt-e_Jac coordinate of Jacobian molecule reference point
// [min,max]_m = range of ipar-wrt-e_Jac molecule m in Jacobian
// [min,max]_fold_m = range of ipar-wrt-e_Jac molecule m which is to folded;
//		      this must be a subrange of [min,max]_m
//
void patch_system::fold_Jacobian(const patch_edge& e_Jac,
				 const patch_edge& e_fold,
				 int iperp,
				 int posn, int min_m, int max_m,
				 int min_fold_m, int max_fold_m,
				 jtutil::array1d<fp>& Jacobian_buffer)
	const
{
// check that [min,max]_fold_m is a subrange of [min,max]_m
assert( min_fold_m >= min_m );
assert( min_fold_m <= max_m );
assert( max_fold_m >= min_m );
assert( max_fold_m <= max_m );

const patch& p = e_fold.my_patch();
assert( e_Jac.my_patch() == p );

const symmetry_ghost_zone& sgz_fold = p.ghost_zone_on_edge(e_fold)
				       .cast_to_symmetry_ghost_zone();

//
// At present we only handle the case show in the comments above,
// where sgz_fold is a local (mirror-image) symmetry, i.e. where
// y is guaranteed to be within the molecule abcde.
//
if (sgz_fold.other_edge() != e_fold)
   then error_exit(ERROR_EXIT,
"***** patch_system::fold_Jacobian()\n"
"        implementation restriction: at present we only handle folding\n"
"        via \"local\" (mirror-image) symmetries!\n"
"        p=\"%s\" e_Jac=\"%s\" e_fold=\"%s\"\n"
"        but sgz_fold.other_edge()=\"%s\" != e_fold\n"
,
		   p.name(), e_Jac.name(), e_fold.name(),
		   sgz_fold.other_edge().name());		/*NOTREACHED*/

	for (int xm = min_fold_m ; xm <= max_fold_m ; ++xm)
	{
	const int x_Jac_ipar = posn + xm;	// x ipar wrt e_Jac
	const int x_fold_iperp = x_Jac_ipar;	// ... == iperp wrt e_fold

	const int y_fold_iperp = sgz_fold.iperp_map_of_iperp(x_fold_iperp);
						// y iperp wrt e_fold
	const int y_Jac_ipar = y_fold_iperp;	// ... == ipar wrt e_Jac
	const int ym = y_Jac_ipar - posn;

	// check that y is indeed within the molecule
	assert( ym >= min_m );
	assert( ym <= max_m );

	// actually "fold" the molecule
	Jacobian_buffer(ym) += Jacobian_buffer(xm);
	}
}

//******************************************************************************

//
// Given that  compute_synchronize_Jacobian()  has been called, this
// function computes a single row of the Jacobian of a given ghost zone,
// *not* taking into account synchronize()'s 3-phase algorithm:
// - It returns the edge to which the y point belongs (the caller can get
//   the patch from this edge).
// - It stores y_iperp and y_posn and min/max ym in the named arguments.
// - It stores the Jacobian elements
//	partial synchronize() gridfn(ghosted_gfn, px, x_iperp, x_ipar)
//	-------------------------------------------------------------
//	     partial gridfn(ghosted_gfn, py, y_iperp, y_posn+ym)
//   in the caller-supplied buffer
//	Jacobian_buffer(ym)
//   for each  ym  in the min/max ym range
//
const patch_edge&
  patch_system::ghost_zone_Jacobian(const ghost_zone& xgz,
				    int x_iperp, int x_ipar,
				    int& y_iperp,
				    int& y_posn, int& min_ym, int& max_ym,
				    jtutil::array1d<fp>& Jacobian_buffer)
	const
{
y_iperp = xgz.Jacobian_y_iperp(x_iperp);

y_posn = xgz.Jacobian_y_ipar_posn(x_iperp, x_ipar);
min_ym = xgz.Jacobian_min_y_ipar_m();
max_ym = xgz.Jacobian_max_y_ipar_m();

	for (int ym = min_ym ; ym <= max_ym ; ++ym)
	{
	Jacobian_buffer(ym) = xgz.Jacobian(x_iperp, x_ipar, ym);
	}

return xgz.other_edge();
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

	  }	// namespace AHFinderDirect
