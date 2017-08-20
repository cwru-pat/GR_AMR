#ifndef AHFD_PATCH_GHOST_ZONE_H
#define AHFD_PATCH_GHOST_ZONE_H
// ghost_zone.hh -- fill in gridfn data in patch ghost zones
// $Header$
//
// ***** design notes for ghost zones *****
// ghost_zone - abstract base class to describe ghost zone of patch
// symmetry_ghost_zone - ... derived class for spacetime-symmetry ghost zone
// interpatch_ghost_zone - ... derived class for interpatch ghost zone
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
//	"../jtutil/cpm_map.hh"
//	"../jtutil/linear_map.hh"
//	"coords.hh"
//	"grid.hh"
//	"fd_grid.hh"
//	"patch.hh"
//	"patch_edge.hh"
//	"patch_interp.hh"
//

// everything in this file is inside this namespace
namespace AHFD
	  {

//*****************************************************************************

//
// ***** design notes for ghost zones *****
//

//
// A  ghost_zone  object describes a patch's ghost zone, and knows how
// to compute gridfns there (we usually speak of "synchronizing" the
// ghost zone or zones) based on either the patch system's symmetry
// or interpolation from a neighboring patch.  ghost_zone is an abstract
// base class, from which we derive two concrete classes:
// * A  symmetry_ghost_zone  object describes a ghost zone which is a
//   (discrete) symmetry of spacetime, either mirror-image or periodic.
//   Such an object knows how to fill in ghost-zone gridfn data from
//   the "other side" of the symmetry.
// * An  interpatch_ghost_zone  object describes a ghost zone which
//   overlaps another patch.  Such an object knows how to get ghost
//   zone gridfn data from the other patch.  More accurately, it gets
//   the data by asking (calling) the appropriate one of the other
//   patch's  patch_interp  objects.
// Every patch has (points to) 4  ghost_zone  objects, one for each of
// the patch's sides.  See the comments in "patch.hh" for a "big picture"
// discussion of patches, patch edges, ghost zones, and patch interpolators.
//

//
// There are some unobvious complications involved in synchronizing
// the ghost zone "corners", i.e. in ghost zone points that are outside
// the nominal grid in *both* coordinates.  There are 3 basic cases here:
// * A corner between two symmetry ghost zones, for example the -x/-y
//   corner in the example below.  In this case it takes *two* sequential
//   symmetry operations to get gridfn data in the corner from the
//   nominal grid.  Symmetry operations commute, so at each point we'll
//   always get the same results independently of in which order we do
//   the symmetry operations.  Computationally, we actually do the operations
//   in both orders, one order's results overwriting the other's, but
//   this doesn't matter (because the results are the same).
// * A corner between two interpatch ghost zones, for example the +x/+y
//   corner in the example below.  In this case we could get the gridfn
//   data by either of two distinct interpolation operations (presumably
//   from two distinct patches), which would in general give slightly
//   different results.  In some ideal world we might do a centered
//   interpolation using data from both patches, but this would be
//   complicated:
//   - it would require a 2-D interpolation
//   - it would require bookkeeping for interpolating from multiple
//     patches within the same ghost zone, indeed for the same ghost
//     zone point
//   At present, we follow a simpler approach: we split the corner down
//   its diagonal,
//	[for the points on the diagonal we make an arbitrary choice;
//	at present this is that they belong to (and get their data via)
//	the rho ghost zone.]
//   and off-center the interpolation as necessary so each ghost-zone
//   point gets data solely from the neighboring patch on its own side.
// * A corner between a symmetry and an interpatch ghost zone, for
//   example the +x/-y or -x/+y corners in the example below.  In this
//   case we first do a symmetry operation in the neighboring patch,
//   then a fully centered interpolation (using the data just obtained
//   from a symmetry operation) to get data in the non-corner part of
//   the interpatch ghost zone.  After the interpatch interpolation,
//   we do a final symmetry operation to get gridfn data in the corner.
//
// In general, then, a ghost zone is rhomboid-shaped: iperp lies in a
// fixed interval, while ipar lies in an interval which may depend on
// iperp.  In general, this shape depends on the type (symmetry vs interpatch)
// of the adjacent ghost zones.
//

//
// To properly handle all the symmetry/interpatch cases described above,
// we use a 3-phase algorithm to synchronize ghost zones:
// Phase 1: Fill in gridfn data at all the non-corner points of symmetry
//	    ghost zones, by using the symmetries to get this data from
//	    its "home patch" nominal grids.
// Phase 2: Fill in gridfn data in all the interpatch ghost zones, by
//	    interpatch interpolating from neighboring patches as described
//	    above.
// Phase 3: Fill in gridfn data at all the corner points of symmetry
//	    ghost zones, by using the symmetries to get this data from
//	    its "home patch" nominal grids or ghost zones.
// Here a given ghost zone corner may be either a full rectangle (so any
// given point is a member of both adjacent corners), or split down its
// diagonal (so any given point is a member of only one corner).  This
// 3-phase algorithm is actually implemented by
//    patch_system::synchronize()
// which in turn calls
//    symmetry_ghost_zone::synchronize()
//    interpatch_ghost_zone::synchronize()
//

//
// For example, consider the +z patch in an octant patch system, with
// the ghost zones being 2 points wide.  The following illustration is
// looking down the z axis, and uses (x,y) for the patch coordinates
// for simplicity:
//
//                    #                                                   //
//                   i+y    i+y    i+y    i+y    i+y    i+y    i+y      //
//   (-2,7) (-1,7)  (0,7)  (1,7)  (2,7)  (3,7)  (4,7)  (5,7)  (6,7)  (7,7)
//    <s-x>  <s-x>    #                                              /i+x
//                    #                                            //
//                   i+y    i+y    i+y    i+y    i+y    i+y      //
//   (-2,6) (-1,6)  (0,6)  (1,6)  (2,6)  (3,6)  (4,6)  (5,6)  (6,6)  (7,6)
//    <s-x>  <s-x>    #                                       /i+x    i+x
//                    #                                     //
//                    #                                   //
//   (-2,5) (-1,5)   2,5)--(1,5)--(2,5)--(3,5)--(4,5)--(5,5)  (6,5)  (7,5)
//     s-x    s-x     #                                  |     i+x    i+x
//                    #                                  |
//                    #                                  |
//   (-2,4) (-1,4)  (0,4)  (1,4)  (2,4)  (3,4)  (4,4)  (5,4)  (6,4)  (7,4)
//     s-x    s-x     #                                  |     i+x    i+x
//                    #                                  |
//                    #                                  |
//   (-2,3) (-1,3)  (0,3)  (1,3)  (2,3)  (3,3)  (4,3)  (5,3)  (6,3)  (7,3)
//     s-x    s-x     #                                  |     i+x    i+x
//                    #                                  |
//                    #                                  |
//   (-2,2) (-1,2)  (0,2)  (1,2)  (2,2)  (3,2)  (4,2)  (5,2)  (6,2)  (7,2)
//     s-x    s-x     #                                  |     i+x    i+x
//                    #                                  |
//                    #                                  |
//   (-2,1) (-1,1)  (0,1)  (1,1)  (2,1)  (3,1)  (4,1)  (5,1)  (6,1)  (7,1)
//     s-x    s-x     #                                  |     i+x    i+x
//                    #                                  |
//                    #                                  |
//  #(-2,0)#(-1,0)##(0,0)##(1,0)##(2,0)##(3,0)##(4,0)##(5,0)##(6,0)##(7,0)
//     s-x    s-x     #                                        i+x    i+x
//                    #
//    <s-y>  <s-y>   s-y    s-y    s-y    s-y    s-y    s-y   <s-y>  <s-y>
//   (-2,-1)(-1,-1) (0,-1) (1,-1) (2,-1) (3,-1) (4,-1) (5,-1) (6,-1) (7,-1)
//    <s-x>  <s-x>    #
//                    #
//    <s-y>  <s-y>   s-y    s-y    s-y    s-y    s-y    s-y   <s-y>  <s-y>
//   (-2,-2)(-1,-2) (0,-2) (1,-2) (2,-2) (3,-2) (4,-2) (5,-2) (6,-2) (7,-2)
//    <s-x>  <s-x>    #
//                    #
//
// For this example,
// * The xz plane and yz plane are marked with ### lines
// * The +z patch's nominal grid is ([0,5],[0,5]), i.e. 0 <= x,y <= 5;
//   its boundary lines are shown with single lines --- and | .
// * The diagonal where we've split corners are marked with // lines.
// * The +z patch's ghost zones are
//	-x: (-1,[-1,7]), (-2,[-2,7]) 
//	+x: (6,[-2,6]), (7,[-2,7])
//	-y: ([-2, 7],[-2,-1])
//	+y: ([-2,5],6), ([-2,6],7)
// * The regions where we will interpolate data from the +z patch are
//	+x: ([ 3,4],[-2,7])
//	+y: ([-2,7],[ 3,4])
//   Note that in both cases the interpolation region includes the points
//   computed by symmetry (in phase 1 of our 3-phase algorithm) on the
//   adjacent edges! There are no interpolation regions inside the -x or
//   -y boundaries, since no interpolation is needed across those boundaries
//   of this patch.
// The diagonal *** line shows the boundary between the +x and +y ghost
// zones.
//
// Our 3-phase algorithm described above thus becomes:
// Phase 1: Fill in gridfn values at points marked with "s-x" below or
//	    "s-y" above via symmetry mirroring across the -x boundary
//	    (yz plane) or -y boundary (xz plane), as described by the
//	    +z patch's -x or -y  symmetry_ghost_zone  object respectively.
// Phase 2: Fill in gridfn values at points marked with "i+x" below or
//	    "i+y" above via interpatch interpolation from the neighboring
//	    patch across the +z patch's +x or +y boundary, as described
//	    by the +z patch's +x or +y  interpatch_ghost_zone  object
//	    respectively.
// Phase 3: Fill in gridfn values at points marked with "<s-x>" below or
//	    "<s-y>" above via symmetry mirroring across the -x boundary
//	    (yz plane) or -y boundary (xz plane), as described by the
//	    +z patch's -x or -y  symmetry_ghost_zone  object respectively.
//

//*****************************************************************************

//
// ghost_zone - abstract base class to describe ghost zone of patch
//
// This is an abstract base class describing a generic patch ghost zone.
// This might represent either of:
// - a discrete symmetry of spacetime (derived class symmetry_ghost_zone)
// - an overlap with another patch (derived class interpatch_ghost_zone)
//

//
// N.b. const qualifiers in ghost_zone and its derived classes refer to
//      the underlying gridfn data.
// 

// forward declarations
class symmetry_ghost_zone;
class interpatch_ghost_zone;
class patch_system;

class	ghost_zone
	{
public:
	//
	// ***** main high-level client interface *****
	//
	// "synchronize" a ghost zone, i.e. update the ghost-zone values
	// of the specified gridfns via the appropriate sequence of
	// symmetry operations and interpatch interpolations
	// (flags specify which part(s) of the ghost zone we want)
	//
	virtual void synchronize(int ghosted_min_gfn, int ghosted_max_gfn,
				 bool want_corners = true,
				 bool want_noncorner = true)
		= 0;


public:
	//
	// ***** Jacobian of synchronize() *****
	//
	// This function computes the Jacobian of the  synchronize()
	// operation into internal buffers; the following functions
	// provide access to that Jacobian.
	//
	// FIXME: should these be moved out into a separate Jacobian
	//        object/class?
	//
	// Note that this function just computes the Jacobian of this
	// ghost zone's  synchronize()  operation -- it does *NOT* take
	// into account the 3-phase synchronization algorithm described
	// in the header comments for this file.  (That's done by
	//  patch_system::synchronize_Jacobian()  and its subfunctions.)
	//
	// n.b. terminology is
	//	partial gridfn at x
	//	-------------------
	//	partial gridfn at y
	//
	virtual void compute_Jacobian(int ghosted_min_gfn, int ghosted_max_gfn,
				      bool want_corners = true,
				      bool want_noncorner = true)
		const
		= 0;

	//
	// The API in the remaining functions implicitly assumes that
	// the Jacobian is independent of  ghosted_gfn , and also that
	// the structure of the Jacobian is such that the set of y points
	// on which a single ghost-zone point depends,
	// - has a single yiperp value (depending on our iperp, of course)
	// - have a contiguous interval of yipar (depending on our iperp
	//   and ipar, of course), whose size is
	//	[or can be taken to be without an unreasonable
	//	amount of zero-padding]
	//   independent of our iperp and ipar; we parameterize this
	//   interval as  yipar = posn+m  where  posn  is determined by
	//   our iperp and ipar, and  m  has a fixed range independent
	//   of our iperp and ipar
	//

	// what is the [min,max] range of m for this ghost zone?
	virtual int Jacobian_min_y_ipar_m() const = 0;
	virtual int Jacobian_max_y_ipar_m() const = 0;

	// what is the iperp of the Jacobian y points in their (y) patch?
	virtual int Jacobian_y_iperp(int x_iperp) const = 0;

	// what is the  posn  value of the y points in this Jacobian row?
	virtual int Jacobian_y_ipar_posn(int x_iperp, int x_ipar) const = 0;

	// what is the Jacobian
	//	partial synchronize() px.gridfn(ghosted_gfn, x_iperp, x_ipar)
	//	-------------------------------------------------------------
	//	   partial py.gridfn(ghosted_gfn, y_iperp, y_posn+y_ipar_m)
	// where
	//	y_iperp = Jacobian_y_iperp(x_iperp)
	//	y_posn = Jacobian_y_ipar_posn(x_iperp, x_ipar)
	virtual fp Jacobian(int x_iperp, int x_ipar, int y_ipar_m) const = 0;


public:
	//
	// ***** low-level client interface *****
	//

	// to which patch/edge do we belong?
	patch&            my_patch() const { return my_patch_; }
	const patch_edge& my_edge()  const { return my_edge_; }

	// from which patch/edge do we get data?
	patch&            other_patch() const { return other_patch_; }
	const patch_edge& other_edge()  const { return other_edge_; }

	// what type of ghost zone are we?
	bool is_interpatch() const { return  is_interpatch_; }
	bool is_symmetry()   const { return !is_interpatch_; }

	// convenience forwarding functions down to patch_edge::
	bool is_min() const { return my_edge().is_min(); }
	bool is_rho() const { return my_edge().is_rho(); }

	// min/max iperp of the ghost zone
	int min_iperp() const
		{
		return my_patch()
		       .minmax_ang_ghost_zone__min_iperp(is_min(), is_rho());
		}
	int max_iperp() const
		{
		return my_patch()
		       .minmax_ang_ghost_zone__max_iperp(is_min(), is_rho());
		}

	// inner/outer iperp of the ghost zone wrt our patch
	int inner_iperp() const { return is_min() ? max_iperp() : min_iperp(); }
	int outer_iperp() const { return is_min() ? min_iperp() : max_iperp(); }

	// extreme min/max ipar that might possibly be part of this ghost zone
	// (derived classes may actually use a subset of this)
	int extreme_min_ipar() const
		{ return my_edge().min_ipar_with_corners(); }
	int extreme_max_ipar() const
		{ return my_edge().max_ipar_with_corners(); }

	// actual min/max ipar in the ghost zone at a particular iperp
	// (may depend on type of the adjacent ghost zones)
	virtual int min_ipar(int iperp) const = 0;
	virtual int max_ipar(int iperp) const = 0;

	// point membership predicate
	bool is_in_ghost_zone(int iperp, int ipar)
		const
		{
		// n.b. don't test ipar until we're sure iperp is in range!
		return (iperp >= min_iperp()) && (iperp <= max_iperp())
		       && (ipar >= min_ipar(iperp))
		       && (ipar <= max_ipar(iperp));
		}

	// adjacent ghost zones to our min/max corners
	const ghost_zone& min_par_adjacent_ghost_zone() const
		{
		return my_patch()
		       .ghost_zone_on_edge( my_edge().min_par_adjacent_edge() );
		}
	const ghost_zone& max_par_adjacent_ghost_zone() const
		{
		return my_patch()
		       .ghost_zone_on_edge( my_edge().max_par_adjacent_edge() );
		}

	//
	// ***** safely cast to derived classes *****
	//

	// assert that gz is of specified type,
	// then static_cast to derive type
	const symmetry_ghost_zone& cast_to_symmetry_ghost_zone() const;
	      symmetry_ghost_zone& cast_to_symmetry_ghost_zone();
	const interpatch_ghost_zone& cast_to_interpatch_ghost_zone() const;
	      interpatch_ghost_zone& cast_to_interpatch_ghost_zone();

	//
	// ***** constructor, finish setup, destructor *****
	//
protected:
	// ... values for  is_interpatch_in  constructor argument
	//     FIXME: these should really be bool, but then we couldn't
	//            use the "enum hack" for in-class constants
	enum {
	     ghost_zone_is_symmetry = false,
	     ghost_zone_is_interpatch = true // no comma
	     };

	// constructor
	// ... only used in implementing our derived classes;
	//     the rest of the world constructs our derived classes instead
	ghost_zone(const patch_edge& my_edge_in,
		   const patch_edge& other_edge_in,
		   bool is_interpatch_in)
		: my_patch_(my_edge_in.my_patch()),
		  my_edge_(my_edge_in),
		  other_patch_(other_edge_in.my_patch()),
		  other_edge_(other_edge_in),
		  is_interpatch_(is_interpatch_in)
		{ }
public:
	// assert() that ghost zone is fully setup:
	// defined here ==> no-op
	// symmetry ghost zone ==> unchanged ==> no-op
	// interpatch ghost zone ==> check consistency of this and the
	//			     other patch's ghost zones and
	//			     patch_interp objects
	virtual void assert_fully_setup() const { }

	// destructor must be virtual to allow destruction
	// of derived classes via ptr/ref to this class
	virtual ~ghost_zone() { }

private:
	// we forbid copying and passing by value
	// by declaring the copy constructor and assignment operator
	// private, but never defining them (either here or in derived classes)
	ghost_zone(const ghost_zone& rhs);
	ghost_zone& operator=(const ghost_zone& rhs);

private:
	patch& my_patch_;
	const patch_edge& my_edge_;
	patch& other_patch_;
	const patch_edge& other_edge_;
	const bool is_interpatch_;
	};

//*****************************************************************************

//
// symmetry_ghost_zone - derived class for spacetime-symmetry ghost zone
//
// In practice, there are two types of spacetime symmetry ghost zone:
// mirror symmetry and periodic symmetry.  However, it turns out that the
// code needed to handle periodic BCs is basically a superset of that
// needed to handle mirror symmetries, so this class represents a generic
// symmetry ghost zone which may be of either type, and once constructed
// doesn't distinguish between the two.
//
// In general, a symmetry ghost zone implies that there's a 1-1 mapping
// between ghost zone points of this patch, and (a subset of the) interior
// points of this or another patch.  If tensors are involved (this isn't
// used at present in the horizon finder), there's also a corresponding
// 1-1 mapping between (angular) tensor components.
//
// A mirror-symmetry ghost zone is specified by (the constructor arguments)
// - a patch edge
// - the (fp) perp coordinate of the mirror plane
// The mapping of ghost zone points is thus "just" the mirror imaging of
// iperp across the symmetry plane within this same patch.  (The mapping
// leaves ipar invariant.)
//
// A periodic-symmetry ghost zone is specified by (the constructor arguments)
// - a patch edge (specifies the ghost zone)
// - the patch edge to which the ghost zone is to be mapped
// - a pair of ipar coordinates, one on this edge and one on the other edge,
//   which map into each other
// - the sign of the ipar mapping (does increasing ipar on this edge map to
//   increasing or decreasing ipar on the other edge?)
// The mapping of ghost zone points is the periodic mapping; this may map
// the ghost zone points to interior points of either this same patch or a
// different one.
//
// In general, the symmetry mapping of ghost zone points is of the form
//	(iperp, ipar) --> (const +/- iperp, const +/- ipar)
// The iperp mapping is always in the direction
//	outside the patch --> inside the patch
// while the ipar mapping might have either sign.
// If there are tensors, the corresponding mapping of tensor components is
//	(index_perp, index_par) --> (+/-) (+/-) (index_perp, index_par)
// (that is, the two +/- signs are multiplied).

//
// Since all the member functions are  const , a  symmetry_ghost_zone
// object is effectively always  const .
//
class	symmetry_ghost_zone
	: public ghost_zone
	{
public:
	//
	// ***** main high-level client interface *****
	//
	// "synchronize" a ghost zone, i.e. update the ghost-zone values
	// of the specified gridfns via the appropriate symmetry operations
	// (flags specify which part(s) of the ghost zone we want)
	//
	void synchronize(int ghosted_min_gfn, int ghosted_max_gfn,
			 bool want_corners = true,
			 bool want_noncorner = true);

	//
	// ***** Jacobian of synchronize() *****
	//
	// n.b. terminology is
	//      partial gridfn at x
	//      -------------------
	//      partial gridfn at y
	//

	// allocate internal buffers, compute Jacobian
	// ... this function is a no-op in this class
	void compute_Jacobian(int ghosted_min_gfn, int ghosted_max_gfn,
			      bool want_corners = true,
			      bool want_noncorner = true)
		const
		{ }

	// what is the [min,max] range of m for this ghost zone?
	int Jacobian_min_y_ipar_m() const { return 0; }
	int Jacobian_max_y_ipar_m() const { return 0; }

	// what is the oiperp of the Jacobian points (= iperp in their patch)?
	virtual int Jacobian_y_iperp(int x_iperp) const
		{ return iperp_map_of_iperp(x_iperp); }

	// what is the  posn  value of the points in this Jacobian row?
	int Jacobian_y_ipar_posn(int x_iperp, int x_ipar) const
		{ return ipar_map_of_ipar(x_ipar); }

	// what is the Jacobian
	//	partial synchronize() px.gridfn(ghosted_gfn, x_iperp, x_ipar)
	//	-------------------------------------------------------------
	//	   partial py.gridfn(ghosted_gfn, y_iperp, y_posn+y_ipar_m)
	// where
	//	y_iperp = Jacobian_y_iperp(x_iperp)
	//	y_posn = Jacobian_y_ipar_posn(x_iperp, x_ipar)
	fp Jacobian(int x_iperp, int x_ipar, int y_ipar_m) const
		{ return (y_ipar_m == 0) ? 1.0 : 0.0; }


	//
	// ***** low-level client interface *****
	//

	// symmetry-map coordinates
	int iperp_map_of_iperp(int iperp) const
		{ return iperp_map_->map(iperp); }
	int ipar_map_of_ipar(int ipar) const
		{ return ipar_map_->map(ipar); }
	fp fp_sign_of_iperp_map() const
		{ return iperp_map_->fp_sign(); }
	fp fp_sign_of_ipar_map() const
		{ return ipar_map_->fp_sign(); }

	// min/max ipar of the ghost zone
	// ... we always include the corners
	//     (cf. the example at the start of this file)
	int min_ipar(int iperp) const { return extreme_min_ipar(); }
	int max_ipar(int iperp) const { return extreme_max_ipar(); }

	//
	// ***** constructors, destructor *****
	//
public:
	// constructor for mirror-symmetry ghost zone
	symmetry_ghost_zone(const patch_edge& my_edge_in);

	// constructor for periodic-symmetry ghost zone
	// ... ipar mapping specified by giving sample point and mapping sign
	symmetry_ghost_zone
	    (const patch_edge& my_edge_in, const patch_edge& other_edge_in,
	     int my_edge_sample_ipar,      int other_edge_sample_ipar,
	     bool ipar_map_is_plus);

	~symmetry_ghost_zone();

private:
        // we forbid copying and passing by value
        // by declaring the copy constructor and assignment operator
        // private, but never defining them
	symmetry_ghost_zone(const symmetry_ghost_zone& rhs);
	symmetry_ghost_zone& operator=(const symmetry_ghost_zone& rhs);

private:
	// symmetry mappings for (iperp,ipar)
	// ... we own these objects
	const jtutil::cpm_map<fp>* iperp_map_;
	const jtutil::cpm_map<fp>* ipar_map_;
	};

//*****************************************************************************

//
// interpatch_ghost_zone - derived class for interpatch ghost zone of a patch
//
// A ghost_zone object maps (my_iperp,my_ipar) coordinates to the other
// patch's (other_iperp,other_par) coordinates, then calls the other patch's
// patch_interp object to interpolate the other patch's data to those
// coordinates.
//
// Note that as described in the "design notes for ghost zones"
// comments above,  interpatch_ghost_zone  objects are constructed in
// the 2nd and 3rd phase of the overall construction process described
// at the comments at the start of "patch.hh"
// [done by our constructor]
// - set up the object itslf and its links to/from the patches and
//   their edges
// [done by  finish_setup()]
// - set up the interpatch mapping information, data pointers, and
//   interpolation result buffer
// - construct the  patch_interp  object to interpolate from the other
//   patch, and save a pointer to it
//

class patch_interp;

class	interpatch_ghost_zone
	: public ghost_zone
	{
public:
	//
	// ***** main high-level client interface *****
	//
	// "synchronize" a ghost zone, i.e. update the ghost-zone
	// values of the specified gridfns via the appropriate
	// interpatch interpolations
	// (flags specify which part(s) of the ghost zone we want)
	//
	// ... the present implementation only supports the case where
	//     both flags are set
	//
	void synchronize(int ghosted_min_gfn, int ghosted_max_gfn,
			 bool want_corners = true,
			 bool want_noncorner = true);

	//
	// ***** Jacobian of synchronize() *****
	//
	// n.b. terminology is
	//      partial gridfn at x
	//      -------------------
	//      partial gridfn at y
	//

	// allocate internal buffers, compute Jacobian
	//
	// ... the present implementation only supports the case where
	//     both flags are set
	//
	void compute_Jacobian(int ghosted_min_gfn, int ghosted_max_gfn,
			      bool want_corners = true,
			      bool want_noncorner = true)
		const;

	// what is the [min,max] range of m for this ghost zone?
	int Jacobian_min_y_ipar_m() const { return Jacobian_min_y_ipar_m_; }
	int Jacobian_max_y_ipar_m() const { return Jacobian_max_y_ipar_m_; }

	// what is the iperp of the Jacobian y points in their (y) patch?
	// ... the ipar row of grid points is actually the same, so
	//     we just have to translate x_iperp to the y patch's coordinates
	int Jacobian_y_iperp(int x_iperp) const { return other_iperp(x_iperp); }

	// what is the  posn  value of the y points in this Jacobian row?
	int Jacobian_y_ipar_posn(int x_iperp, int x_ipar) const
		{
		assert(Jacobian_y_ipar_posn_ != NULL);
		const int y_iperp = Jacobian_y_iperp(x_iperp);
		return (*Jacobian_y_ipar_posn_)(y_iperp, x_ipar);
		}

	// what is the Jacobian
	//	partial synchronize() px.gridfn(ghosted_gfn, x_iperp, x_ipar)
	//	-------------------------------------------------------------
	//	   partial py.gridfn(ghosted_gfn, y_iperp, y_posn+y_ipar_m)
	// where
	//	y_iperp = Jacobian_y_iperp(x_iperp)
	//	y_posn = Jacobian_y_ipar_posn(x_iperp, x_ipar)
	fp Jacobian(int x_iperp, int x_ipar, int y_ipar_m) const
		{
		assert(Jacobian_buffer_ != NULL);
		assert(y_ipar_m >= Jacobian_min_y_ipar_m_);
		assert(y_ipar_m <= Jacobian_max_y_ipar_m_);
		const int y_iperp = Jacobian_y_iperp(x_iperp);
		return (*Jacobian_buffer_)(y_iperp, x_ipar, y_ipar_m);
		}


	//
	// ***** low-level client interface *****
	//

public:
	// check consistency of this and the other patch's ghost zones
	// and patch_interp objects
	void assert_fully_setup() const;

	// min/max ipar of the ghost zone for specified iperp
	// with possibly "triangular" corners depending on the type
	// (symmetry vs interpatch) of the adjacent ghost zones
	// (cf. comments & example at the start of this file)
	int min_ipar(int iperp) const;
	int max_ipar(int iperp) const;

	// convert our iperp --> other patch's iperp
	int other_iperp(int iperp) const
		{
		assert(other_iperp_ != NULL);
		return other_iperp_->map(iperp);
		}

	//
	// ***** constructor, finish setup, destructor *****
	//
public:
	interpatch_ghost_zone(const patch_edge& my_edge_in,
			      const patch_edge& other_edge_in,
			      int patch_overlap_width);

	// finish setup (requires adjacent-side ghost_zone objects
	// to exist, though not to have finish_setup() called):
	// - setup par coordinate mapping information
	// - setup interpatch interpolator data pointers & result buffer
	// - create patch_interp object to interpolate from *other* patch
	void finish_setup(int interp_handle, int interp_par_table_handle);

	~interpatch_ghost_zone();

private:
        // we forbid copying and passing by value
        // by declaring the copy constructor and assignment operator
        // private, but never defining them
	interpatch_ghost_zone(const interpatch_ghost_zone& rhs);
	interpatch_ghost_zone& operator=(const interpatch_ghost_zone& rhs);

private:
	//
	// all the remaining pointers are initialized to NULL pointers
	// in our constructor, then finally allocated and set up by
	// finish_setup() or compute_Jacobian() as appropriate
	//
	// FIXME: should these be moved out into a separate object/class
	//        for the interp stuff and/or another one for the Jacobian?
	//

	// see comment in "patch_interp.hh" for why this is "const"
	const patch_interp* other_patch_interp_;

	// other patch's iperp coordinates of our ghost zone points
	// ... maps my_iperp --> other_iperp
	jtutil::cpm_map<fp>* other_iperp_;

	// min/max values of other patch's iperp coordinates
	// of our ghost zone points
	int min_other_iperp_, max_other_iperp_;

	// [min,max]_ipar used at each other_iperp
	// ... we will pass these arrays by reference
	//     to the other patch's patch_interp object
	// ... index is (other_iperp)
	jtutil::array1d<int>* min_ipar_used_;
	jtutil::array1d<int>* max_ipar_used_;

	// other patch's (fp) parallel coordinates of our ghost zone points
	// ... we will pass this array by reference
	//     to the other patch's patch_interp object
	//     using my_ipar as the patch_interp's parindex
	// ... subscripts are (other_iperp, my_ipar)
	jtutil::array2d<fp>* other_par_;

	// buffer into which the other patch's patch_interp object
	// will store the interpolated gridfn values
	// ... we will pass this array by reference
	//     to the other patch's patch_interp object
	//     using my_ipar as the patch_interp's parindex
	// ... subscripts are (gfn, other_iperp,my_ipar)
	jtutil::array3d<fp>* interp_result_buffer_;

	//
	// stuff computed by  compute_Jacobian()
	//
	// n.b. terminology is
	//      partial gridfn at x
	//      -------------------
	//      partial gridfn at y
	//
        mutable int Jacobian_min_y_ipar_m_;
	mutable int Jacobian_max_y_ipar_m_;

	// other patch's y ipar posn for a Jacobian row
	// ... subscripts are (oiperp, ipar)
	mutable jtutil::array2d<CCTK_INT>* Jacobian_y_ipar_posn_;

	// Jacobian values
	// ... subscripts are (y_iperp, x_ipar, y_ipar_m)
	mutable jtutil::array3d<fp>* Jacobian_buffer_;
	};

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif
