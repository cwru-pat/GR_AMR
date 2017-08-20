// ghost_zone.cc -- fill in gridfn data in patch ghost zones
// $Header$

//
// ghost_zone::cast_to_symmetry_ghost_zone
// ghost_zone::cast_to_interpatch_ghost_zone
//
// symmetry_ghost_zone::symmetry_ghost_zone (mirror symmetry)
// symmetry_ghost_zone::symmetry_ghost_zone (periodic BC)
// symmetry_ghost_zone::~symmetry_ghost_zone
// symmetry_ghost_zone::synchronize
//
// interpatch_ghost_zone::interpatch_ghost_zone
// interpatch_ghost_zone::~interpatch_ghost_zone
// interpatch_ghost_zone::[min,max]_ipar
// interpatch_ghost_zone::finish_setup
// interpatch_ghost_zone::assert_fully_setup
// interpatch_ghost_zone::synchronize
// interpatch_ghost_zone::compute_Jacobian
//

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

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
namespace AHFinderDirect
	  {

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// These functions verify (assert()) that a ghost zone is indeed of
// the specified type, then static_cast to the appropriate derived class.
//

const symmetry_ghost_zone& ghost_zone::cast_to_symmetry_ghost_zone()
	const
{
assert( is_symmetry() );
return static_cast<const symmetry_ghost_zone &>(*this);
}

symmetry_ghost_zone& ghost_zone::cast_to_symmetry_ghost_zone()
{
assert( is_symmetry() );
return static_cast<symmetry_ghost_zone &>(*this);
}

//**************************************

const interpatch_ghost_zone& ghost_zone::cast_to_interpatch_ghost_zone()
	const
{
assert( is_interpatch() );
return static_cast<const interpatch_ghost_zone &>(*this);
}

interpatch_ghost_zone& ghost_zone::cast_to_interpatch_ghost_zone()
{
assert( is_interpatch() );
return static_cast<interpatch_ghost_zone &>(*this);
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function constructs a mirror-symmetry ghost zone object
//
symmetry_ghost_zone::symmetry_ghost_zone(const patch_edge& my_edge_in)
	: ghost_zone(my_edge_in,
		     my_edge_in,	// other edge == my edge
		     ghost_zone_is_symmetry)
{
// iperp_map: i --> (i of ghost zone) - i
iperp_map_ = new jtutil::cpm_map<fp>(min_iperp(), max_iperp(),
				     my_edge_in.fp_grid_outer_iperp());

// ipar_map_: identity map
ipar_map_ = new jtutil::cpm_map<fp>(extreme_min_ipar(), extreme_max_ipar());
}

//******************************************************************************

//
// This function constructs a periodic-symmetry ghost zone object.
//
symmetry_ghost_zone::symmetry_ghost_zone
	(const patch_edge& my_edge_in, const patch_edge& other_edge_in,
	 int my_edge_sample_ipar,      int other_edge_sample_ipar,
	 bool ipar_map_is_plus)
	: ghost_zone(my_edge_in,
		     other_edge_in,
		     ghost_zone_is_symmetry)
{
//
// perpendicular map
//
const fp fp_my_period_plane_iperp    = my_edge()   .fp_grid_outer_iperp();
const fp fp_other_period_plane_iperp = other_edge().fp_grid_outer_iperp();

// iperp mapping must be outside --> inside
// i.e. if both edges have iperp as the same min/max "direction",
//	then the mapping is  iperp increasing --> iperp decreasing
//      (i.e. the map's sign is -1)
const bool is_iperp_map_plus
	= ! (my_edge().is_min() == other_edge().is_min());
iperp_map_ = new jtutil::cpm_map<fp>(min_iperp(), max_iperp(),
				     fp_my_period_plane_iperp,
				     fp_other_period_plane_iperp,
				     is_iperp_map_plus);

//
// parallel map
//
ipar_map_ = new jtutil::cpm_map<fp>(extreme_min_ipar(), extreme_max_ipar(),
				    my_edge_sample_ipar, other_edge_sample_ipar,
				    ipar_map_is_plus);
}

//******************************************************************************

//
// This function destroys a  symmetry_ghost_zone  object.
//
symmetry_ghost_zone::~symmetry_ghost_zone()
{
delete ipar_map_;
delete iperp_map_;
}

//******************************************************************************

//
// This function "synchronizes" a ghost zone, i.e. it updates the
// ghost-zone values of the specified gridfns via the appropriate
// symmetry operations.The flags specify which part(s) of the ghost zone
// we want.
//
void symmetry_ghost_zone::synchronize(int ghosted_min_gfn, int ghosted_max_gfn,
				      bool want_corners /* = true */,
				      bool want_noncorner /* = true */)
{
	for (int gfn = ghosted_min_gfn ; gfn <= ghosted_max_gfn ; ++gfn)
	{
	for (int iperp = min_iperp() ; iperp <= max_iperp() ; ++iperp)
	{
	for (int ipar = min_ipar(iperp) ; ipar <= max_ipar(iperp) ; ++ipar)
	{
	// do we want to do this point?
	if (! my_edge().ipar_is_in_selected_part(want_corners, want_noncorner,
						 ipar) )
	   then continue;				// *** LOOP CONTROL ***

	const int sym_iperp = iperp_map_of_iperp(iperp);
	const int sym_ipar  = ipar_map_of_ipar  (ipar );
	const int sym_irho = other_edge()
			     .irho_of_iperp_ipar  (sym_iperp,sym_ipar);
	const int sym_isigma = other_edge()
			       .isigma_of_iperp_ipar(sym_iperp,sym_ipar);
	const fp sym_gridfn = other_patch()
			      .ghosted_gridfn(gfn, sym_irho,sym_isigma);

	const int irho   = my_edge().  irho_of_iperp_ipar(iperp,ipar);
	const int isigma = my_edge().isigma_of_iperp_ipar(iperp,ipar);
	my_patch().ghosted_gridfn(gfn, irho,isigma) = sym_gridfn;
	}
	}
	}
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function constructs an  interpatch_ghost_zone  object.
//
interpatch_ghost_zone::interpatch_ghost_zone(const patch_edge& my_edge_in,
					     const patch_edge& other_edge_in,
					     int patch_overlap_width)
	: ghost_zone(my_edge_in,
		     other_edge_in,
		     ghost_zone_is_interpatch),
	  // remaining pointers are all set up properly by finish_setup()
	  other_patch_interp_(NULL),
	  other_iperp_(NULL),
          min_ipar_used_(NULL),	max_ipar_used_(NULL),
	  other_par_(NULL),
	  interp_result_buffer_(NULL),
	  Jacobian_y_ipar_posn_(NULL), Jacobian_buffer_(NULL) // no comma
{
//
// verify that we have the expected relationships between
// this and the other patch's (mu,nu,phi) coordinates:
//

// perp coordinate is common to us and the other patch, so
// ghost zone must be min in one patch, max in the other
if (my_edge().is_min() == other_edge().is_min())
   then error_exit(ERROR_EXIT,
"***** interpatch_ghost_zone::interpatch_ghost_zone:\n"
"        my_patch().name()=\"%s\" my_edge().name()=%s\n"
"        other_patch().name()=\"%s\" other_edge().name()=%s\n"
"        ghost zone must be min in one patch, max in the other!\n"
,
		   my_patch().name(), my_edge().name(),
		   other_patch().name(), other_edge().name());	/*NOTREACHED*/

// coord in common between the two patches must be perp coord in both patches
// and this patch's tau coordinate must be other edge's parallel coordinate
const local_coords::coords_set common_coords_set
	= local_coords::coords_set_not(my_patch().coords_set_rho_sigma()
				       ^
				       other_patch().coords_set_rho_sigma());
if (! (    (common_coords_set == my_edge().coords_set_perp())
	&& (common_coords_set == other_edge().coords_set_perp())
	&& (my_patch().coords_set_tau() == other_edge().coords_set_par())    ) )
   then error_exit(PANIC_EXIT,
"***** interpatch_ghost_zone::interpatch_ghost_zone:\n"
"        (rho,sigma,tau) coordinates don't match up properly\n"
"        between this patch/edge and the other patch/edge!\n"
"        my_patch().name()=\"%s\" my_edge().name()=%s\n"
"        other_patch().name()=\"%s\" other_edge().name()=%s\n"
"        my_patch().coords_set_{rho,sigma,tau}={%s,%s,%s}\n"
"        my_edge().coords_set_{perp,par}={%s,%s}\n"
"        other_patch().coords_set_{rho,sigma,tau}={%s,%s,%s}\n"
"        other_edge().coords_set_{perp,par}={%s,%s}\n"
,
	my_patch().name(), my_edge().name(),
	other_patch().name(), other_edge().name(),
	local_coords::name_of_coords_set(my_patch().coords_set_rho()),
	local_coords::name_of_coords_set(my_patch().coords_set_sigma()),
	local_coords::name_of_coords_set(my_patch().coords_set_tau()),
	local_coords::name_of_coords_set(my_edge().coords_set_perp()),
	local_coords::name_of_coords_set(my_edge().coords_set_par()),
	local_coords::name_of_coords_set(other_patch().coords_set_rho()),
	local_coords::name_of_coords_set(other_patch().coords_set_sigma()),
	local_coords::name_of_coords_set(other_patch().coords_set_tau()),
	local_coords::name_of_coords_set(other_edge().coords_set_perp()),
	local_coords::name_of_coords_set(other_edge().coords_set_par()));
								/*NOTREACHED*/

// perp coordinate must match (mod 2*pi) across the two patches
// after taking into account any overlap
// ... eg patch_overlap_width = 3 would be
//	p   p   p   p   p
//		q   q   q   q   q
//     so the overlap would be (patch_overlap_width-1) * delta
const fp other_overlap
	= (patch_overlap_width-1) * other_edge().perp_map().delta_fp();
const fp other_outer_perp_minus_overlap	// move back inwards into other patch
					// by overlap distance, to get a value
					// that should match our own
					// grid_outer_perp() value
	= other_edge().grid_outer_perp()
	  + (other_edge().is_min() ? + other_overlap : - other_overlap);
if (! local_coords::fuzzy_EQ_ang(my_edge().grid_outer_perp(),
				 other_outer_perp_minus_overlap))
   then error_exit(ERROR_EXIT,
"***** interpatch_ghost_zone::interpatch_ghost_zone:\n"
"        my_patch().name()=\"%s\" my_edge().name()=%s\n"
"        other_patch().name()=\"%s\" other_edge().name()=%s\n"
"        perp coordinate doesn't match (mod 2*pi) across the two patches!\n"
"        my_edge().grid_outer_perp()=%g   <--(compare this)\n"
"        patch_overlap_width=%d other_overlap=%g\n"
"        other_edge.grid_outer_perp()=%g\n"
"        other_outer_perp_minus_overlap=%g   <--(against this)\n"
,
		   my_patch().name(), my_edge().name(),
		   other_patch().name(), other_edge().name(),
		   double(my_edge().grid_outer_perp()),
		   patch_overlap_width, double(other_overlap),
		   double(other_edge().grid_outer_perp()),
		   double(other_outer_perp_minus_overlap));	/*NOTREACHED*/


//
// set up the iperp interpatch coordinate mapping
// (gives other patch's iperp coordinate for interpolation)
//

// compute the iperp --> other_iperp mapping for a sample point;
// ... if the ghost zone is empty, then the sample point will necessarily
//     be out-of-range in the ghost zone, so we use the *unchecked*
//     conversions to avoid errors in this case
// ... we do the computation using the fact that  perp  is the same
//     coordinate in both patches (modulo 2*pi radians = 360 degrees)
const int sample_iperp = outer_iperp();
const fp sample_perp = my_edge().perp_map()
				.fp_of_int_unchecked(sample_iperp);
						// unchecked conversion here!
const fp other_sample_perp = other_patch()
			     .modulo_reduce_ang(other_edge().perp_is_rho(),
						sample_perp);
const fp fp_other_sample_iperp = other_edge()
				 .fp_iperp_of_perp(other_sample_perp);

// verify that this is fuzzily a grid point
if (! jtutil::fuzzy<fp>::is_integer(fp_other_sample_iperp))
   then error_exit(ERROR_EXIT,
"***** interpatch_ghost_zone::interpatch_ghost_zone:\n"
"        my_patch().name()=\"%s\" my_edge().name()=%s\n"
"        other_patch().name()=\"%s\" other_edge().name()=%s\n"
"        sample_iperp=%d sample_perp=%g\n"
"        other_sample_perp=%g fp_other_sample_iperp=%g\n"
"        ==> fp_other_sample_iperp isn't fuzzily an integer!\n"
"        ==> patches aren't commensurate in the perpendicular coordinate!\n"
,
		   my_patch().name(), my_edge().name(),
		   other_patch().name(), other_edge().name(),
		   sample_iperp, double(sample_perp),
		   double(other_sample_perp),
		   double(fp_other_sample_iperp));		/*NOTREACHED*/
const int other_sample_iperp
	= jtutil::round<fp>::to_integer(fp_other_sample_iperp);

// compute the +/- sign (direction) of the iperp --> other_iperp mapping
//
// Since perp is the same in both patches (mod 2*pi radians = 360 degrees),
// the overall +/- sign is just the product of the signs of the two individual
// iperp <--> perp mappings.
//
// ... signs encoded as (floating-point) +/- 1.0
const double iperp_map_sign_pm1
	=   jtutil::signum(    my_edge().perp_map().delta_fp() )
	  * jtutil::signum( other_edge().perp_map().delta_fp() );
// ... signs encoded as is_plus bool flag
const bool is_iperp_map_plus = (iperp_map_sign_pm1 > 0.0);

// now we finally know enough to set up the other_iperp(iperp)
// coordinate mapping
other_iperp_ = new jtutil::cpm_map<fp>(min_iperp(), max_iperp(),
				       sample_iperp, other_sample_iperp,
				       is_iperp_map_plus);
}

//******************************************************************************

//
// this function destroys an  interpatch_ghost_zone  object.
//
interpatch_ghost_zone::~interpatch_ghost_zone()
{
delete Jacobian_buffer_;
delete Jacobian_y_ipar_posn_;
delete interp_result_buffer_;
delete other_par_;
delete max_ipar_used_;
delete min_ipar_used_;
delete other_iperp_;
delete other_patch_interp_;
}

//******************************************************************************

//
// These functions compute the [min,max] ipar of the ghost zone for
// a given iperp, taking into account how we treat the corners
// (cf. the example in the header comments in "ghost_zone.hh"):
//
// If an adjacent ghost zone is symmetry,
//    we do not include that corner;
// If an adjacent ghost zone is interpatch,
//    we include up to the diagonal line, and if we are a rho ghost zone,
//    then also the diagonal line itself.  E.g. For the example in the
//    header comments "ghost_zone.hh", the +x ghost zone includes (6,6),
//    (7,6), and (7,7), while the +y ghost zone includes (6,7)
//
// ... in the following 2 functions,
//     the  iabs()  term includes the diagonal,
//     so we must remove the diagonal for !is_rho,
//     i.e. add 1 to min_ipar and subtract 1 from max_ipar
//
int interpatch_ghost_zone::min_ipar(int iperp) const
{
return min_par_adjacent_ghost_zone().is_symmetry()
       ? my_edge().min_ipar_without_corners()
       : my_edge().min_ipar_without_corners()
	 - iabs(iperp - my_edge().nominal_grid_outer_iperp())
	 + (is_rho() ? 0 : 1);
}

int interpatch_ghost_zone::max_ipar(int iperp) const
{
return max_par_adjacent_ghost_zone().is_symmetry()
       ? my_edge().max_ipar_without_corners()
       : my_edge().max_ipar_without_corners()
	 + iabs(iperp - my_edge().nominal_grid_outer_iperp())
	 - (is_rho() ? 0 : 1);
}

//******************************************************************************

//
// This function finishes the construction/setup of an  interpatch_ghost_zone
// object.  It
// - sets up the par coordinate mapping information
// - sets up the interpatch interpolator data pointer and result arrays
// - constructs the patch_interp object to interpolate from the *other* patch
//
// We use our ipar as the patch_interp's parindex.
//
void interpatch_ghost_zone::finish_setup(int interp_handle,
					 int interp_par_table_handle)
{
min_other_iperp_ = jtutil::min(other_iperp(min_iperp()),
			       other_iperp(max_iperp()));
max_other_iperp_ = jtutil::max(other_iperp(min_iperp()),
			       other_iperp(max_iperp()));


//
// set up arrays giving actual [min,max] ipar that we'll use
// at each other_iperp (later on we will pass these arrays to the
// other patch's  patch_interp  object, with ipar being parindex there
//
min_ipar_used_ = new jtutil::array1d<int>(min_other_iperp_, max_other_iperp_);
max_ipar_used_ = new jtutil::array1d<int>(min_other_iperp_, max_other_iperp_);
	  {
	for (int iperp = min_iperp() ; iperp <= max_iperp() ; ++iperp)
	{
	(*min_ipar_used_)(other_iperp(iperp)) = min_ipar(iperp);
	(*max_ipar_used_)(other_iperp(iperp)) = max_ipar(iperp);
	}
	  }


//
// set up array giving other patch's par coordinate for interpolation
//

other_par_ = new jtutil::array2d<fp>(min_other_iperp_, max_other_iperp_,
				     extreme_min_ipar(), extreme_max_ipar());

	  {
	for (int iperp = min_iperp() ; iperp <= max_iperp() ; ++iperp)
	{
	for (int ipar = min_ipar(iperp); ipar <= max_ipar(iperp) ; ++ipar)
	{
	// compute the  other_par corresponding to  (iperp,ipar)
	// ... here we use the fact (which we verified in our constructor)
	//     that other edge's parallel coordinate == our tau coordinate
	//     (at least modulo 2*pi radians = 360 degrees)
	const fp perp  = my_edge().perp_of_iperp(iperp);
	const fp par   = my_edge().par_of_ipar(ipar);

	const fp rho   = my_edge().  rho_of_perp_par(perp, par);
	const fp sigma = my_edge().sigma_of_perp_par(perp, par);

	const fp tau   = my_patch().tau_of_rho_sigma(rho, sigma);
	const fp other_par = other_patch()
			     .modulo_reduce_ang(other_edge().par_is_rho(), tau);

	(*other_par_)(other_iperp(iperp),ipar) = other_par;
	}
	}
	  }


//
// set up interpolation result buffer
//
interp_result_buffer_
	= new jtutil::array3d<fp>(my_patch().ghosted_min_gfn(),
					  my_patch().ghosted_max_gfn(),
				  min_other_iperp_, max_other_iperp_,
				  extreme_min_ipar(), extreme_max_ipar());

//
// construct the patch_interp object to interpolate from the *other* patch
// ... the patch_interp should use gridfn data from it's (the other patch's)
//     min/max par ghost zones if those (adjacent) adjacent ghost zones
//     are symmetry, but not if they're interpatch,
//     cf the header comments in "ghost_zone.hh"
//
const ghost_zone& other_ghost_zone = other_patch()
				     .ghost_zone_on_edge(other_edge());
const bool ok_to_use_min_par_ghost_zone
	= other_ghost_zone.min_par_adjacent_ghost_zone()
			  .is_symmetry()
	  ? true
	  : false;
const bool ok_to_use_max_par_ghost_zone
	= other_ghost_zone.max_par_adjacent_ghost_zone()
			  .is_symmetry()
	  ? true
	  : false;
other_patch_interp_ = new patch_interp(other_edge(),
				       min_other_iperp_, max_other_iperp_,
				       *min_ipar_used_, *max_ipar_used_,
				       *other_par_,
				       ok_to_use_min_par_ghost_zone,
				       ok_to_use_max_par_ghost_zone,
				       interp_handle, interp_par_table_handle);
}


//******************************************************************************

//
// This function asserts() that
// - we have a patch_interp object
// - our and the patch_interp object's notions of the "other patch" agree
// - the other patch has an interpatch ghost zone on this edge
// - the other patch's interpatch ghost zone on this edge,
//   points back to our patch
//
void interpatch_ghost_zone::assert_fully_setup() const
{
assert(other_patch_interp_ != NULL);
assert(other_patch() == other_patch_interp_->my_patch());
assert( other_patch()
	.ghost_zone_on_edge(other_edge())
	.is_interpatch() );
assert( my_patch() == other_patch()
		      .ghost_zone_on_edge(other_edge())
		      .other_patch() );
}

//******************************************************************************

//
// This function "synchronizes" a ghost zone, i.e. it updates the
// ghost-zone values of the specified gridfns via the appropriate
// interpatch interpolations.
//
// The flags specify which part(s) of the ghost zone we want, but
// the present implementation only supports the case where all the
// flags are  true , i.e. we want the entire ghost zone.
//
void interpatch_ghost_zone::synchronize
	(int ghosted_min_gfn, int ghosted_max_gfn,
	 bool want_corners /* = true */,
	 bool want_noncorner /* = true */)
{
// make sure the caller wants the entire ghost zone
if (! (want_corners && want_noncorner))
   then error_exit(ERROR_EXIT,
"***** interpatch_ghost_zone::synchronize():\n"
"        we only support operating on the *entire* ghost zone,\n"
"        but we were passed flags specifying a proper subset!\n"
"        want_corners=(int)%d want_noncorner=(int)%d\n"
,
		   want_corners, want_noncorner);		/*NOTREACHED*/

// do the interpolation into our result buffer
other_patch_interp_->interpolate(ghosted_min_gfn, ghosted_max_gfn,
				 *interp_result_buffer_);

// store the results back into our gridfns
    for (int gfn = ghosted_min_gfn ; gfn <= ghosted_max_gfn ; ++gfn)
    {
	for (int iperp = min_iperp() ; iperp <= max_iperp() ; ++iperp)
	{
	const int oiperp = other_iperp(iperp);

	for (int ipar = min_ipar(iperp) ; ipar <= max_ipar(iperp) ; ++ipar)
	{
	int irho   = my_edge().  irho_of_iperp_ipar(iperp,ipar);
	int isigma = my_edge().isigma_of_iperp_ipar(iperp,ipar);
	my_patch().ghosted_gridfn(gfn, irho,isigma)
		= (*interp_result_buffer_)(gfn, oiperp,ipar);
	}
	}
    }
}

//******************************************************************************

//
// This function allocates the internal buffers for the Jacobian, and
// computes that Jacobian
//	    partial synchronize gridfn(ghosted_gfn, iperp, ipar)
//	------------------------------------------------------------
//	partial other patch gridfn(ghosted_gfn, oiperp, posn+ipar_m)
// where
//	oiperp = Jacobian_oiperp(iperp)
//	posn = Jacobian_oipar_posn(iperp, ipar)
// into the internal buffers.
//
void interpatch_ghost_zone::compute_Jacobian
	(int ghosted_min_gfn, int ghosted_max_gfn,
	 bool want_corners /* = true */,
	 bool want_noncorner /* = true */)
	const
{
// make sure the caller wants the entire ghost zone
if (! (want_corners && want_noncorner))
   then error_exit(ERROR_EXIT,
"***** interpatch_ghost_zone::compute_Jacobian():\n"
"        we only support operating on the *entire* ghost zone,\n"
"        but we were passed flags specifying a proper subset!\n"
"        want_corners=(int)%d want_noncorner=(int)%d\n"
,
		   want_corners, want_noncorner);		/*NOTREACHED*/

assert(other_patch_interp_ != NULL);
other_patch_interp_->verify_Jacobian_sparsity_pattern_ok();

other_patch_interp_->molecule_minmax_ipar_m(Jacobian_min_y_ipar_m_,
					    Jacobian_max_y_ipar_m_);

if (Jacobian_y_ipar_posn_ == NULL)
   then Jacobian_y_ipar_posn_ = new jtutil::array2d<CCTK_INT>
				     (min_other_iperp_, max_other_iperp_,
				      extreme_min_ipar(), extreme_max_ipar());
other_patch_interp_->molecule_posn(*Jacobian_y_ipar_posn_);

if (Jacobian_buffer_ == NULL)
   then Jacobian_buffer_
		= new jtutil::array3d<fp>
			 (min_other_iperp_, max_other_iperp_,
			  extreme_min_ipar(), extreme_max_ipar(),
			  Jacobian_min_y_ipar_m_, Jacobian_max_y_ipar_m_);
other_patch_interp_->Jacobian(*Jacobian_buffer_);
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

	  }	// namespace AHFinderDirect
