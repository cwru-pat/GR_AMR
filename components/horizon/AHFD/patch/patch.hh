#ifndef AHFD_PATCH_PATCH_H
#define AHFD_PATCH_PATCH_H
// patch.hh -- describes a coordinate/grid patch
// $Header$
//
// ***** how patch boundaries are handled *****
// patch - abstract base class to describe a coordinate/grid patch
//
// z_patch - derived class for a +/- z patch
// x_patch - derived class for a +/- x patch
// y_patch - derived class for a +/- y patch
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
//

// everything in this file is inside this namespace
namespace AHFinderDirect
	  {

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

//
// ***** how patch boundaries are handled *****
//

//
// Basically, we handle patch boundaries using the usual "ghost zone"
// technique, interpolating values from neighboring patches as necessary.
//
// In more detail, we use the following interrelated types of objects
// to handle patch boundaries:
//
// A  patch_edge  object represents the basic geometry of a min/max
// rho/sigma side of a patch, i.e. it provides which-side-am-I predicates,
// coordinate conversions between (perp,par) and (rho,sigma), etc.
// Every patch has (points to) 4  patch_edge  objects, one for each of
// the patch's sides.
//
// A  ghost_zone  object describes a patch's ghost zone, and knows how
// to fill in gridfns there based on either the patch system's symmetry
// or interpolation from a neighboring patch.  ghost_zone is an abstract
// base class, from which we derive two classes:
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
// the patch's sides.
//
// A  patch_interp  object does the actual interpolation of data from
// within a patch (for filling in data in another patch's ghost zone).
// A  patch_interp  object points to the patch and patch_edge where it
// will be interpolating.
//
// For example, suppose we have two patches p and q with a common
// angular boundary.  Then the desired network of pointers looks like
// this (omitting the  patch_edge  objects for simplicity):
//
// +-----+                                                       +-----+
// |     | <--> p.interpatch_ghost_zone ---> q.patch_interp ---> |     |
// |  p  |                                                       |  q  |
// |     | <--- p.patch_interp <--- q.interpatch_ghost_zone <--> |     |
// +-----+                                                       +-----+
//
// Because of the mutual pointers, we can't easily construct (say) p's
// interpatch_ghost_zone until after q itself has been constructed, and
// vice versa.  Moreover, the  patch_interp::  constructor needs the
// adjacent-side  ghost_zone  objects to already exist, and it needs to
// know the iperp range of the interpolation region, which can only be
// computed from the adjacent-patch  interpatch_ghost_zone  object.
//
// The solution adopted here is to use a 3-phase algorithm, ultimately
// driven by the patch_system constructor:
// * The patch constructors themselves construct the  patch_edge  objects
//   and links them to/from the patches.
// * The patch_system constructor calls the appropriate functions
//	patch::create_mirror_symmetry_ghost_zone()
//	patch::create_periodic_symmetry_ghost_zone()
//	patch::create_interpatch_ghost_zone()
//   to construct the  ghost_zone  objects and link them to/from the
//   patches.
// * The patch_system constructor calls the functions
//	interpatch_ghost_zone::finish_setup()
//   to finish setting up the  interpatch_ghost_zone  objects, construct
//   the other patch's  patch_interp  objects, and finish linking the
//    interpatch_ghost_zone  objects to the  patch_interp  objects.
//

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

//
// patch - abstract base class to describe a generic coordinate/grid patch
//

//
// There are 3 types of patches, z, x, and y.  Each type uses two of
// (mu,nu,phi) as its angular coordinates (rho,sigma); the remaining
// "unused" one of (mu,nu,phi) is tau.
//
//	z patch ==> (rho,sigma) = (mu,nu)    tau = phi
//	x patch ==> (rho,sigma) = (nu,phi)   tau = mu
//	y patch ==> (rho,sigma) = (mu,phi)   tau = nu
//

// forward declarations
class patch_edge;
class ghost_zone;
class symmetry_ghost_zone;
class interpatch_ghost_zone;
class patch_interp;
class patch_system;

//
// const qualifiers refer to the gridfn values
//
class	patch
	: public fd_grid
	{
	//
	// ***** patch system, type, and coordinate metadata *****
	//
public:

	// to which patch system do we belong?
	patch_system& my_patch_system() const
		{ return my_patch_system_; }

	// each patch has a unique 0-origin small-integer patch number,
	// usually denoted  pn
	int patch_number() const { return patch_number_; }

	// each patch has a unique human-readable patch name for debugging etc
	const char* name() const { return name_; }	// typically "+z" etc

	// are we a +[xyz] or -[xyz] patch?
	bool is_plus() const { return is_plus_; }

	// ... values for the  is_plus_in  constructor argument
	//     FIXME: these should really be bool, but then we couldn't
	//            use the "enum hack" for in-class constants
	enum { patch_is_plus = true, patch_is_minus = false };

	// are we a (+/-) x or y or z patch?
	// ... n.b. type is `char' because this is handy for both
	//	    switch() and human-readable printing
	char ctype() const { return ctype_; }		// 'z' or 'x' or 'y'

	// are two patches really the same patch?
	// n.b. this does *not* compare any of the gridfn data!
	bool operator==(const patch& other_patch) const
		{ return this == &other_patch; }
	bool operator!=(const patch& other_patch) const
		{ return ! operator==(other_patch); }

	// (rho,sigma,tau) coordinates as singleton coordinate sets
	local_coords::coords_set coords_set_rho() const
		{ return coords_set_rho_; } 
	local_coords::coords_set coords_set_sigma() const
		{ return coords_set_sigma_; }
	local_coords::coords_set coords_set_tau() const
		{ return coords_set_tau_; }

	// {rho,sigma} coordinate set
	local_coords::coords_set coords_set_rho_sigma() const
		{ return coords_set_rho() | coords_set_sigma(); }

	// (rho,sigma) coordinates as human-readable character strings
	// (for labelling output files etc)
	virtual const char* name_of_rho() const = 0;
	virtual const char* name_of_sigma() const = 0;


	//
	// ***** (rho,sigma,tau) coordinates *****
	//
public:

	// convert (rho,sigma) --> tau
	virtual fp tau_of_rho_sigma(fp rho, fp sigma) const = 0;

	// convert (rho,sigma) --> (mu,nu,phi)
	virtual fp  mu_of_rho_sigma(fp rho, fp sigma) const = 0;
	virtual fp  nu_of_rho_sigma(fp rho, fp sigma) const = 0;
	virtual fp phi_of_rho_sigma(fp rho, fp sigma) const = 0;

	// convert (rho,sigma) <--> usual polar spherical (theta,phi)
	virtual void theta_phi_of_rho_sigma(fp rho, fp sigma,
					    fp& ps_theta, fp& ps_phi)
		const = 0;
	virtual void rho_sigma_of_theta_phi(fp ps_theta, fp ps_phi,
					    fp& rho, fp& sigma)
		const = 0;

	// convert (r,rho,sigma) <--> local (x,y,z)
	virtual void xyz_of_r_rho_sigma(fp r, fp rho, fp sigma,
					fp& x, fp& y, fp& z)
		const = 0;
	virtual fp   rho_of_xyz(fp x, fp y, fp z) const = 0;
	virtual fp sigma_of_xyz(fp x, fp y, fp z) const = 0;

	// convert (rho,sigma) --> direction cosines (xcos,ycos,zcos)
	//                         with respect to the local coordinate system
	virtual void xyzcos_of_rho_sigma(fp rho, fp sigma,
					 fp& xcos, fp& ycos, fp& zcos)
		const = 0;

	// partial (x,y,z) / partial (rho,sigma)
	virtual void partial_xyz_wrt_r_rho_sigma
	   (fp r, fp rho, fp sigma,
	    fp& partial_x_wrt_r, fp& partial_x_wrt_rho, fp& partial_x_wrt_sigma,
	    fp& partial_y_wrt_r, fp& partial_y_wrt_rho, fp& partial_y_wrt_sigma,
	    fp& partial_z_wrt_r, fp& partial_z_wrt_rho, fp& partial_z_wrt_sigma)
		const = 0;

	// partial (rho,sigma) / partial (x,y,z)
	virtual fp partial_rho_wrt_x(fp x, fp y, fp z) const = 0;
	virtual fp partial_rho_wrt_y(fp x, fp y, fp z) const = 0;
	virtual fp partial_rho_wrt_z(fp x, fp y, fp z) const = 0;
	virtual fp partial_sigma_wrt_x(fp x, fp y, fp z) const = 0;
	virtual fp partial_sigma_wrt_y(fp x, fp y, fp z) const = 0;
	virtual fp partial_sigma_wrt_z(fp x, fp y, fp z) const = 0;

	// partial^2 (rho,sigma) / partial (xx,xy,xz,yy,yz)
	virtual fp partial2_rho_wrt_xx(fp x, fp y, fp z) const = 0;
	virtual fp partial2_rho_wrt_xy(fp x, fp y, fp z) const = 0;
	virtual fp partial2_rho_wrt_xz(fp x, fp y, fp z) const = 0;
	virtual fp partial2_rho_wrt_yy(fp x, fp y, fp z) const = 0;
	virtual fp partial2_rho_wrt_yz(fp x, fp y, fp z) const = 0;
	virtual fp partial2_rho_wrt_zz(fp x, fp y, fp z) const = 0;
	virtual fp partial2_sigma_wrt_xx(fp x, fp y, fp z) const = 0;
	virtual fp partial2_sigma_wrt_xy(fp x, fp y, fp z) const = 0;
	virtual fp partial2_sigma_wrt_xz(fp x, fp y, fp z) const = 0;
	virtual fp partial2_sigma_wrt_yy(fp x, fp y, fp z) const = 0;
	virtual fp partial2_sigma_wrt_yz(fp x, fp y, fp z) const = 0;
	virtual fp partial2_sigma_wrt_zz(fp x, fp y, fp z) const = 0;

	// compute (rho,sigma) 2-D induced metric from 3-D xyz metric
	// as per p.33 of my apparent horizon finding notes
	// ... returns Jacobian of (rho,sigma) 2-D induced metric
	fp rho_sigma_metric(fp r, fp rho, fp sigma,
			    fp partial_surface_r_wrt_rho,
			    fp partial_surface_r_wrt_sigma,
			    fp g_xx, fp g_xy, fp g_xz,
				     fp g_yy, fp g_yz,
					      fp g_zz,
			    fp& g_rho_rho, fp& g_rho_sigma,
					   fp& g_sigma_sigma)
		const;

	// plotting coordinates (dpx,dpy)
	// ... character string describing how (dpx,dpy) are
	//     defined in terms of (mu,nu,phi), eg "90 - drho = 90 - dphi"
	//     (used for labelling output files)
	virtual const char* name_of_dpx() const = 0;
	virtual const char* name_of_dpy() const = 0;
	// ... (irho,isimga) --> (px,py)
	virtual fp dpx_of_rho_sigma(fp rho, fp sigma) const = 0;
	virtual fp dpy_of_rho_sigma(fp rho, fp sigma) const = 0;


	//
	// ***** line/surface integrals *****
	//
public:

	//
	// The following enum describes the integration methods supported
	// by  integrate_gridfn() .
	//
	// For convenience of exposition we describe the methods as if for
	// 1-D integration, but  integrate_gridfn()  actually does 2-D
	// (surface) integration over the patch.
	//
	// Suppose we're computing $\int_{x_0}^{x^N} f(x) \, dx$, using the
	// equally spaced integration points $f_0$, $f_1$, \dots, $f_N$,
	// spaced $\Delta x$ apart.  Then the integration methods are as
	// follows, with the convention that $\langle X \rangle$ denotes
	// indefinite repetition of the "X" terms, depending on N:
	//
	enum	integration_method
		{
		// Trapezoid rule
		// ... character-string name "trapezoid" or "trapezoid rule"
		// ... 2nd order accurate for smooth functions
		// ... requires N >= 1
		// $$
		// \Delta x \left[
		//            \half f_0
		//          + \langle
		//            f_k
		//            \rangle
		//          + \half f_N
		//          \right]
		// $$
		integration_method__trapezoid,

		// Simpson's rule
		// ... character-string name "Simpson" or "Simpson's rule"
		// ... 4th order accurate for smooth functions
		// ... requires N >= 2 and N even
		// $$
		// \Delta x \left[
		//            \frac{1}{3} f_0
		//          + \frac{4}{3} f_1
		//          + \langle
		//            \frac{2}{3} f_{2k} + \frac{4}{3} f_{2k+1}
		//            \rangle
		//          + \frac{1}{3} f_N
		//          \right]
		// $$
		integration_method__Simpson,

		// Simpson's rule, variant form
		// ... characgter-string name "Simpson (variant)"
		//     or "Simpson's rule (variant)"
		// ... described in Numerical Recipes 1st edition (4.1.14)
		// ... 4th order accurate for smooth functions
		// ... requires N >= 7
		// $$
		// \Delta x \left[
		//            \frac{17}{48} f_0
		//          + \frac{59}{48} f_1
		//          + \frac{43}{48} f_2
		//          + \frac{49}{48} f_3
		//          + \langle
		//            f_k
		//            \rangle
		//          + \frac{49}{48} f_{N-3}
		//          + \frac{43}{48} f_{N-2}
		//          + \frac{59}{48} f_{N-1}
		//          + \frac{17}{48} f_N
		//          \right]
		// $$
		integration_method__Simpson_variant,

		// automatic choice of the "best" one of the above methods:
		// ... i.e. choose Simpson's rule or variant if applicable,
		//     otherwise trapezoid rule
		// N == 2	Simpson's rule
		// N == 3	trapezoid rule
		// N == 4	Simpson's rule
		// N == 5	trapezoid rule
		// N == 6	Simpson's rule
		// N >= 7	Simpson's rule, variant form
		integration_method__automatic_choice // no comma here!
		};

	// decode character string name into internal enum
	static
	  enum integration_method
	    decode_integration_method(const char method_string[]);

	// compute the arc length of a surface in the specified plane
	// (must be one of "xy", "xz", or "yz") over the patch's nominal bounds
	// ... error_exit() if  plane  is invalid and/or
	//     the patch doesn't contain that coordinate plane
	virtual fp plane_arc_length(const char plane[],
				    int ghosted_radius_gfn,
				    int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
						  int g_yy_gfn, int g_yz_gfn,
								int g_zz_gfn,
				    enum integration_method method)
		const = 0;

	// ... along the rho direction (i.e. in a dsigma=constant plane
	//     where dsigma is a multiple of 90 degrees)
	fp rho_arc_length(int ghosted_radius_gfn,
			  int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
					int g_yy_gfn, int g_yz_gfn,
						      int g_zz_gfn,
			  enum integration_method method)
		const;
	// ... along the sigma direction (i.e. in a drho=constant plane
	//     where drho is a multiple of 90 degrees)
	fp sigma_arc_length(int ghosted_radius_gfn,
			    int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
					  int g_yy_gfn, int g_yz_gfn,
							int g_zz_gfn,
			    enum integration_method method)
		const;

	// compute the surface integral of a gridfn over the patch's
	// nominal area,
	//	$\int f(\rho,\sigma) \, dA$
	//		= \int f(\rho,\sigma) \sqrt{|J|} \, d\rho \, d\sigma$
	// where $J$ is the Jacobian of $(x,y,z)$ with respect to $(rho,sigma)
	// ... integration method selected by  method  argument
	// ... src gridfn may be either nominal-grid or ghosted-grid
	//     (n.b. in the latter case the integral is still done
	//           only over the patch's nominal area)
	fp integrate_gridfn(int unknown_src_gfn,
			    int ghosted_radius_gfn,
			    int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
					  int g_yy_gfn, int g_yz_gfn,
							int g_zz_gfn,
			    enum integration_method method)
		const;

	fp integrate_gridpoint(int unknown_src_gfn,
                               int ghosted_radius_gfn,
                               int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
                                             int g_yy_gfn, int g_yz_gfn,
                                                           int g_zz_gfn,
                               enum integration_method method,
                               int irho, int isigma)
		const;

	// compute integration coefficient $c_i$ where
	// $\int_{x_0}^{x_N} f(x) \, dx
	//	\approx \Delta x \, \sum_{i=0}^N c_i f(x_i)$
private:
	static
	  fp integration_coeff(enum integration_method method, int N, int i);


	//
	// ***** patch edges ****
	//
public:
	const patch_edge& min_rho_patch_edge() const
		{ return min_rho_patch_edge_; }
	const patch_edge& max_rho_patch_edge() const
		{ return max_rho_patch_edge_; }
	const patch_edge& min_sigma_patch_edge() const
		{ return min_sigma_patch_edge_; }
	const patch_edge& max_sigma_patch_edge() const
		{ return max_sigma_patch_edge_; }
	const patch_edge& minmax_ang_patch_edge(bool want_min, bool want_rho)
		const
		{
		return want_min ? (want_rho ? min_rho_patch_edge()
					    : min_sigma_patch_edge())
				: (want_rho ? max_rho_patch_edge()
					    : max_sigma_patch_edge());
		}

	// find which patch edge is adjacent to neighboring patch q,
	// or error_exit() if it's not actually a neighboring patch
	// ... computation done using only (rho,sigma) coordinate sets
	//     and min/max dang bounds ==> ok to use in setting up ghost zones
	// ... patch_overlap_width = number of grid points (grid spacings
	//     in the perpendicular direction) these patches' nominal grids
	//     overlap,
	//     ... if this is nonzero, then these patches must have
	//         the *same* grid spacing in the perpendicular direction
	//     ... e.g. delta_dang = 5, this patch max_dang = 50,
	//         other patch min_dang = 40 ==> patch_overlap_width = 3
	//		p   p   p   p   p
	//			q   q   q   q   q
	const patch_edge& edge_adjacent_to_patch(const patch& q,
						 int patch_overlap_width = 0)
		const;


	//
	// ***** ghost zones *****
	//
public:
	ghost_zone& min_rho_ghost_zone() const
		{
		assert(min_rho_ghost_zone_ != NULL);
		return *min_rho_ghost_zone_;
		}
	ghost_zone& max_rho_ghost_zone() const
		{
		assert(max_rho_ghost_zone_ != NULL);
		return *max_rho_ghost_zone_;
		}
	ghost_zone& min_sigma_ghost_zone() const
		{
		assert(min_sigma_ghost_zone_ != NULL);
		return *min_sigma_ghost_zone_;
		}
	ghost_zone& max_sigma_ghost_zone() const
		{
		assert(max_sigma_ghost_zone_ != NULL);
		return *max_sigma_ghost_zone_;
		}
	ghost_zone& minmax_rho_ghost_zone(bool want_min)
		const
		{
		return want_min ? min_rho_ghost_zone()
				: max_rho_ghost_zone();
		}
	ghost_zone& minmax_sigma_ghost_zone(bool want_min)
		const
		{
		return want_min ? min_sigma_ghost_zone()
				: max_sigma_ghost_zone();
		}

	ghost_zone& minmax_ang_ghost_zone(bool want_min, bool want_rho)
		const
		{
		return want_rho ? minmax_rho_ghost_zone(want_min)
				: minmax_sigma_ghost_zone(want_min);
		}

	ghost_zone& ghost_zone_on_edge(const patch_edge &e) const;

	// which of the two ghost zones at a specified corner,
	// contains a specified point?
	ghost_zone& corner_ghost_zone_containing_point
		(bool rho_is_min, bool sigma_is_min,	// specifies corner
		 int irho, int isigma)			// specifies point
		const;

	// which ghost zone contains a specified noncorner point?
	ghost_zone& ghost_zone_containing_noncorner_point(int irho, int isigma)
		const;


	//
	// ***** set up ghost zones
	//
public:

	// assert() that this ghost zone hasn't been set up yet,
	// then set it up as mirror-symmetry
	void create_mirror_symmetry_ghost_zone(const patch_edge& edge);

	// assert() that this ghost zone hasn't been set up yet,
	// then set it up as periodic-symmetry
	void create_periodic_symmetry_ghost_zone
		(const patch_edge& my_edge, const patch_edge& other_edge,
		 bool ipar_map_is_plus);

	// assert() that this ghost zone hasn't been set up yet,
	// then set it up as interpatch
	// ... this only sets up ghost zone in skeletal form; use
	//     interpatch_ghost_zone::finish_setup()  to complete
	//     the setup process
	void create_interpatch_ghost_zone
		(const patch_edge& my_edge, const patch_edge& other_edge,
		 int patch_overlap_width);

	// assert() that all ghost zones
	// are fully setup
	void assert_all_ghost_zones_fully_setup() const;

private:
	// helper function for setup_*_ghost_zone():
	// assert() that ghost zone pointer on specified edge is NULL
	// (i.e. that we haven't already setup this ghost zone),
	// then assign new value to it
	void set_ghost_zone(const patch_edge& edge, ghost_zone* gzp);


	//
	// ***** constructor, destructor, et al *****
	//
protected:
	// ... used only from derived classes
	// ... doesn't set up ghost zone info, since this depends on
	//     knowing our neighbouring patches, which might not exist yet
	// ... saves a pointer to name_in[], so this should have a
	//     lifetime at least as long as that of this object
	patch(patch_system &my_patch_system_in, int patch_number_in,
	      const char name_in[], bool is_plus_in, char ctype_in,
	      local_coords::coords_set coords_set_rho_in,
	      local_coords::coords_set coords_set_sigma_in,
	      local_coords::coords_set coords_set_tau_in,
	      const grid_arrays::grid_array_pars& grid_array_pars_in,
	      const grid::grid_pars& grid_pars_in);
public:
	// destructor must be virtual to allow destruction
	// of derived classes via ptr/ref to this class
	virtual ~patch();

private:
	// we forbid copying and passing by value
	// by declaring the copy constructor and assignment operator
	// private, but never defining them
	patch(const patch& rhs);
	patch& operator=(const patch& rhs);


	//
	// ***** data members *****
	//
private:

	// type/coordinate metadata
	patch_system &my_patch_system_;
	const int patch_number_;
	const char* name_;
	const bool is_plus_;
	const char ctype_;
	const local_coords::coords_set coords_set_rho_,
				       coords_set_sigma_,
				       coords_set_tau_;

	// edges
	const patch_edge& min_rho_patch_edge_;
	const patch_edge& max_rho_patch_edge_;
	const patch_edge& min_sigma_patch_edge_;
	const patch_edge& max_sigma_patch_edge_;

	// ghost zones
	// ... pointers are set to NULL by ctor,
	//     reset to non-NULL by set_ghost_zone(), which is called by
	//	  create_mirror_symmetry_ghost_zone()
	//	  create_periodic_symmetry_ghost_zone()
	//	  create_interpatch_ghost_zone()
	ghost_zone* min_rho_ghost_zone_;
	ghost_zone* max_rho_ghost_zone_;
	ghost_zone* min_sigma_ghost_zone_;
	ghost_zone* max_sigma_ghost_zone_;
	};

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

//
// This class describes a +/- z patch.  It doesn't define any new
// functions not already present in  class patch ; it "just" defines
// non-virtual versions of all the pure virtual functions defined there.
//
//	z patch ==> (rho,sigma) = (mu,nu)    tau = phi
//
class	z_patch
	: public patch
	{
public:
	// human-readable names of (rho,sigma)
	const char* name_of_rho() const { return "mu"; }
	const char* name_of_sigma() const { return "nu"; }

	// convert (rho,sigma) --> tau
	fp tau_of_rho_sigma(fp rho, fp sigma) const
		{ return local_coords::phi_of_mu_nu(rho,sigma); }

	// convert (rho,sigma) --> (mu,nu,phi)
	fp mu_of_rho_sigma(fp rho, fp sigma) const { return rho; }
	fp nu_of_rho_sigma(fp rho, fp sigma) const { return sigma; }
	fp phi_of_rho_sigma(fp rho, fp sigma) const
		{ return local_coords::phi_of_mu_nu(rho,sigma); }

	// convert (rho,sigma) <--> usual polar spherical (theta,phi)
	void theta_phi_of_rho_sigma(fp rho, fp sigma, fp& ps_theta, fp& ps_phi)
		const
		{
		local_coords::theta_phi_of_mu_nu(rho,sigma, ps_theta,ps_phi);
		}
	void rho_sigma_of_theta_phi(fp ps_theta, fp ps_phi, fp& rho, fp& sigma)
		const
		{
		local_coords::mu_nu_of_theta_phi(ps_theta,ps_phi, rho,sigma);
		}

	// convert (r,rho,sigma) <--> (x,y,z)
	void xyz_of_r_rho_sigma(fp r, fp rho, fp sigma, fp& x, fp& y, fp& z)
		const
		{ local_coords::xyz_of_r_mu_nu(r,rho,sigma, x,y,z); }
	fp rho_of_xyz(fp x, fp y, fp z) const
		{ return modulo_reduce_rho(local_coords::mu_of_yz(y,z)); }
	fp sigma_of_xyz(fp x, fp y, fp z) const
		{ return modulo_reduce_sigma(local_coords::nu_of_xz(x,z)); }

	// convert (rho,sigma) --> direction cosines (xcos,ycos,zcos)
	//                         with respect to the local coordinate system
	void xyzcos_of_rho_sigma(fp rho, fp sigma,
				 fp& xcos, fp& ycos, fp& zcos)
		const
		{ local_coords::xyzcos_of_mu_nu(rho,sigma, xcos,ycos,zcos); }

	// partial (x,y,z) / partial (rho,sigma)
	void partial_xyz_wrt_r_rho_sigma
	   (fp r, fp rho, fp sigma,
	    fp& partial_x_wrt_r, fp& partial_x_wrt_rho, fp& partial_x_wrt_sigma,
	    fp& partial_y_wrt_r, fp& partial_y_wrt_rho, fp& partial_y_wrt_sigma,
	    fp& partial_z_wrt_r, fp& partial_z_wrt_rho, fp& partial_z_wrt_sigma)
		const
		{
		local_coords::partial_xyz_wrt_r_mu_nu
		   (r, rho, sigma,
		    partial_x_wrt_r, partial_x_wrt_rho, partial_x_wrt_sigma,
		    partial_y_wrt_r, partial_y_wrt_rho, partial_y_wrt_sigma,
		    partial_z_wrt_r, partial_z_wrt_rho, partial_z_wrt_sigma);
		}

	// partial (rho,sigma) / partial (x,y,z)
	fp partial_rho_wrt_x(fp x, fp y, fp z) const { return 0.0; }
	fp partial_rho_wrt_y(fp x, fp y, fp z) const
		{ return local_coords::partial_mu_wrt_y(y,z); }
	fp partial_rho_wrt_z(fp x, fp y, fp z) const
		{ return local_coords::partial_mu_wrt_z(y,z); }
	fp partial_sigma_wrt_x(fp x, fp y, fp z) const
		{ return local_coords::partial_nu_wrt_x(x,z); }
	fp partial_sigma_wrt_y(fp x, fp y, fp z) const { return 0.0; }
	fp partial_sigma_wrt_z(fp x, fp y, fp z) const
		{ return local_coords::partial_nu_wrt_z(x,z); }

	// partial^2 (rho,sigma) / partial (xx,xy,xz,yy,yz)
	fp partial2_rho_wrt_xx(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_rho_wrt_xy(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_rho_wrt_xz(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_rho_wrt_yy(fp x, fp y, fp z) const
		{ return local_coords::partial2_mu_wrt_yy(y,z); }
	fp partial2_rho_wrt_yz(fp x, fp y, fp z) const
		{ return local_coords::partial2_mu_wrt_yz(y,z); }
	fp partial2_rho_wrt_zz(fp x, fp y, fp z) const
		{ return local_coords::partial2_mu_wrt_zz(y,z); }
	fp partial2_sigma_wrt_xx(fp x, fp y, fp z) const
		{ return local_coords::partial2_nu_wrt_xx(x,z); }
	fp partial2_sigma_wrt_xy(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_sigma_wrt_xz(fp x, fp y, fp z) const
		{ return local_coords::partial2_nu_wrt_xz(x,z); }
	fp partial2_sigma_wrt_yy(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_sigma_wrt_yz(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_sigma_wrt_zz(fp x, fp y, fp z) const
		{ return local_coords::partial2_nu_wrt_zz(x,z); }

	// plotting coordinates (px,py)
	// ... character string describing how (dpx,dpy) are
	//     defined in terms of (mu,nu,phi), eg "90 - drho = 90 - dphi"
	//     (used for labelling output files)
	const char* name_of_dpx() const
		{ return "dsigma = dnu"; }
	const char* name_of_dpy() const
		{ return is_plus() ? "drho = dmu" : "180 - drho = 180 - dmu"; }
	// ... (irho,isimga) --> (px,py)
	fp dpx_of_rho_sigma(fp rho, fp sigma) const
		{ return jtutil::degrees_of_radians(sigma); }
	fp dpy_of_rho_sigma(fp rho, fp sigma) const
		{
		const fp drho = jtutil::degrees_of_radians(rho);
		return is_plus() ? drho : 180.0 - drho;
		}

	// compute the arc length of a surface in the specified plane
	// (must be one of "xz" or "yz") over the patch's nominal bounds
	// ... error_exit() if  plane  is invalid
	fp plane_arc_length(const char plane[],
			    int ghosted_radius_gfn,
			    int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
					  int g_yy_gfn, int g_yz_gfn,
							int g_zz_gfn,
			    enum integration_method method)
		const;

	// constructor, destructor
	z_patch(patch_system &my_patch_system_in, int patch_number_in,
		const char* name_in, bool is_plus_in,
		const grid_arrays::grid_array_pars& grid_array_pars_in,
		const grid::grid_pars& grid_pars_in);
	~z_patch() { }

private:
        // we forbid copying and passing by value
        // by declaring the copy constructor and assignment operator
        // private, but never defining them
	z_patch(const z_patch& rhs);
	z_patch& operator=(const z_patch& rhs);
	};

//*****************************************************************************

//
// This class describes a +/- x patch.  It doesn't define any new
// functions not already present in  class patch ; it "just" defines
// non-virtual versions of all the pure virtual functions defined there.
//
//	x patch ==> (rho,sigma) = (nu,phi)   tau = mu
//
class	x_patch
	: public patch
	{
public:
	// human-readable names of (rho,sigma)
	const char* name_of_rho() const { return "nu"; }
	const char* name_of_sigma() const { return "phi"; }

	// convert (rho,sigma) --> tau
	fp tau_of_rho_sigma(fp rho, fp sigma) const
		{ return local_coords::mu_of_nu_phi(rho, sigma); }

	// convert (rho,sigma) --> (mu,nu,phi)
	fp nu_of_rho_sigma(fp rho, fp sigma) const { return rho; }
	fp phi_of_rho_sigma(fp rho, fp sigma) const { return sigma; }
	fp mu_of_rho_sigma(fp rho, fp sigma) const
		{ return local_coords::mu_of_nu_phi(rho, sigma); }

	// convert (rho,sigma) <--> usual polar spherical (theta,phi)
	void theta_phi_of_rho_sigma(fp rho, fp sigma, fp& ps_theta, fp& ps_phi)
		const
		{
		local_coords::theta_phi_of_nu_phi(rho, sigma, ps_theta, ps_phi);
		}
	void rho_sigma_of_theta_phi(fp ps_theta, fp ps_phi, fp& rho, fp& sigma)
		const
		{
		local_coords::nu_phi_of_theta_phi(ps_theta, ps_phi, rho, sigma);
		}

	// convert (r,rho,sigma) <--> (x,y,z)
	void xyz_of_r_rho_sigma(fp r, fp rho, fp sigma, fp& x, fp& y, fp& z)
		const
		{ local_coords::xyz_of_r_nu_phi(r, rho, sigma, x, y, z); }
	fp rho_of_xyz(fp x, fp y, fp z) const
		{ return modulo_reduce_rho(local_coords::nu_of_xz(x, z)); }
	fp sigma_of_xyz(fp x, fp y, fp z) const
		{ return modulo_reduce_sigma(local_coords::phi_of_xy(x, y)); }

	// convert (rho,sigma) --> direction cosines (xcos,ycos,zcos)
	//                         with respect to the local coordinate system
	void xyzcos_of_rho_sigma(fp rho, fp sigma,
				 fp& xcos, fp& ycos, fp& zcos)
		const
		{ local_coords::xyzcos_of_nu_phi(rho,sigma, xcos,ycos,zcos); }

	// partial (x,y,z) / partial (rho,sigma)
	void partial_xyz_wrt_r_rho_sigma
	   (fp r, fp rho, fp sigma,
	    fp& partial_x_wrt_r, fp& partial_x_wrt_rho, fp& partial_x_wrt_sigma,
	    fp& partial_y_wrt_r, fp& partial_y_wrt_rho, fp& partial_y_wrt_sigma,
	    fp& partial_z_wrt_r, fp& partial_z_wrt_rho, fp& partial_z_wrt_sigma)
		const
		{
		local_coords::partial_xyz_wrt_r_nu_phi
		   (r, rho, sigma,
		    partial_x_wrt_r, partial_x_wrt_rho, partial_x_wrt_sigma,
		    partial_y_wrt_r, partial_y_wrt_rho, partial_y_wrt_sigma,
		    partial_z_wrt_r, partial_z_wrt_rho, partial_z_wrt_sigma);
		}

	// partial (rho,sigma) / partial (x,y,z)
	fp partial_rho_wrt_x(fp x, fp y, fp z) const
		{ return local_coords::partial_nu_wrt_x(x,z); }
	fp partial_rho_wrt_y(fp x, fp y, fp z) const { return 0.0; }
	fp partial_rho_wrt_z(fp x, fp y, fp z) const
		{ return local_coords::partial_nu_wrt_z(x,z); }
	fp partial_sigma_wrt_x(fp x, fp y, fp z) const
		{ return local_coords::partial_phi_wrt_x(x,y); }
	fp partial_sigma_wrt_y(fp x, fp y, fp z) const
		{ return local_coords::partial_phi_wrt_y(x,y); }
	fp partial_sigma_wrt_z(fp x, fp y, fp z) const { return 0.0; }

	// partial^2 (rho,sigma) / partial (xx,xy,xz,yy,yz)
	fp partial2_rho_wrt_xx(fp x, fp y, fp z) const
		{ return local_coords::partial2_nu_wrt_xx(x,z); }
	fp partial2_rho_wrt_xy(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_rho_wrt_xz(fp x, fp y, fp z) const
		{ return local_coords::partial2_nu_wrt_xz(x,z); }
	fp partial2_rho_wrt_yy(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_rho_wrt_yz(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_rho_wrt_zz(fp x, fp y, fp z) const
		{ return local_coords::partial2_nu_wrt_zz(x,z); }
	fp partial2_sigma_wrt_xx(fp x, fp y, fp z) const
		{ return local_coords::partial2_phi_wrt_xx(x,y); }
	fp partial2_sigma_wrt_xy(fp x, fp y, fp z) const
		{ return local_coords::partial2_phi_wrt_xy(x,y); }
	fp partial2_sigma_wrt_xz(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_sigma_wrt_yy(fp x, fp y, fp z) const
		{ return local_coords::partial2_phi_wrt_yy(x,y); }
	fp partial2_sigma_wrt_yz(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_sigma_wrt_zz(fp x, fp y, fp z) const { return 0.0; }

	// plotting coordinates (px,py)
	// ... character string describing how (dpx,dpy) are
	//     defined in terms of (mu,nu,phi), eg "90 - drho = 90 - dphi"
	//     (used for labelling output files)
	const char* name_of_dpx() const { return "drho = dnu"; }
	const char* name_of_dpy() const
		{
		return is_plus() ? "dsigma = dphi"
				 : "180 - dsigma = 180 - dphi";
		}
	// ... (irho,isimga) --> (px,py)
	fp dpx_of_rho_sigma(fp rho, fp sigma) const
		{ return jtutil::degrees_of_radians(rho); }
	fp dpy_of_rho_sigma(fp rho, fp sigma) const
		{
		const fp dsigma = jtutil::degrees_of_radians(sigma);
		return is_plus() ? dsigma : 180.0 - dsigma;
		}

	// compute the arc length of a surface in the specified plane
	// (must be one of "xy" or "xz") over the patch's nominal bounds
	// ... error_exit() if  plane  is invalid
	fp plane_arc_length(const char plane[],
			    int ghosted_radius_gfn,
			    int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
					  int g_yy_gfn, int g_yz_gfn,
							int g_zz_gfn,
			    enum integration_method method)
		const;

	// constructor, destructor
	x_patch(patch_system &my_patch_system_in, int patch_number_in,
		const char* name_in, bool is_plus_in,
		const grid_arrays::grid_array_pars& grid_array_pars_in,
		const grid::grid_pars& grid_pars_in);
	~x_patch() { }

private:
        // we forbid copying and passing by value
        // by declaring the copy constructor and assignment operator
        // private, but never defining them
	x_patch(const x_patch& rhs);
	x_patch& operator=(const x_patch& rhs);
	};

//*****************************************************************************

//
// This class describes a +/- y patch.  It doesn't define any new
// functions not already present in  class patch ; it "just" defines
// non-virtual versions of all the pure virtual functions defined there.
//
//	y patch ==> (rho,sigma) = (mu,phi)   tau = nu
//
class	y_patch
	: public patch
	{
public:
	// human-readable names of (rho,sigma)
	const char* name_of_rho() const { return "mu"; }
	const char* name_of_sigma() const { return "phi"; }

	// convert (rho,sigma) --> tau
	fp tau_of_rho_sigma(fp rho, fp sigma) const
		{ return local_coords::nu_of_mu_phi(rho, sigma); }

	// convert (rho,sigma) --> (mu,nu,phi)
	fp mu_of_rho_sigma(fp rho, fp sigma) const { return rho; }
	fp phi_of_rho_sigma(fp rho, fp sigma) const { return sigma; }
	fp nu_of_rho_sigma(fp rho, fp sigma) const
		{ return local_coords::nu_of_mu_phi(rho, sigma); }

	// convert (rho,sigma) <--> usual polar spherical (theta,phi)
	void theta_phi_of_rho_sigma(fp rho, fp sigma, fp& ps_theta, fp& ps_phi)
		const
		{
		local_coords::theta_phi_of_mu_phi(rho, sigma, ps_theta, ps_phi);
		}
	void rho_sigma_of_theta_phi(fp ps_theta, fp ps_phi, fp& rho, fp& sigma)
		const
		{
		local_coords::mu_phi_of_theta_phi(ps_theta, ps_phi, rho, sigma);
		}

	// convert (r,rho,sigma) <--> (x,y,z)
	void xyz_of_r_rho_sigma(fp r, fp rho, fp sigma, fp& x, fp& y, fp& z)
		const
		{ local_coords::xyz_of_r_mu_phi(r, rho, sigma, x, y, z); }
	fp rho_of_xyz(fp x, fp y, fp z) const
		{ return modulo_reduce_rho(local_coords::mu_of_yz(y, z)); }
	fp sigma_of_xyz(fp x, fp y, fp z) const
		{ return modulo_reduce_sigma(local_coords::phi_of_xy(x, y)); }

	// convert (rho,sigma) --> direction cosines (xcos,ycos,zcos)
	//                         with respect to the local coordinate system
	void xyzcos_of_rho_sigma(fp rho, fp sigma,
				 fp& xcos, fp& ycos, fp& zcos)
		const
		{ local_coords::xyzcos_of_mu_phi(rho,sigma, xcos,ycos,zcos); }

	// partial (x,y,z) / partial (rho,sigma)
	void partial_xyz_wrt_r_rho_sigma
	   (fp r, fp rho, fp sigma,
	    fp& partial_x_wrt_r, fp& partial_x_wrt_rho, fp& partial_x_wrt_sigma,
	    fp& partial_y_wrt_r, fp& partial_y_wrt_rho, fp& partial_y_wrt_sigma,
	    fp& partial_z_wrt_r, fp& partial_z_wrt_rho, fp& partial_z_wrt_sigma)
		const
		{
		local_coords::partial_xyz_wrt_r_mu_phi
		   (r, rho, sigma,
		    partial_x_wrt_r, partial_x_wrt_rho, partial_x_wrt_sigma,
		    partial_y_wrt_r, partial_y_wrt_rho, partial_y_wrt_sigma,
		    partial_z_wrt_r, partial_z_wrt_rho, partial_z_wrt_sigma);
		}

	// partial (rho,sigma) / partial (x,y,z)
	fp partial_rho_wrt_x(fp x, fp y, fp z) const { return 0.0; }
	fp partial_rho_wrt_y(fp x, fp y, fp z) const
		{ return local_coords::partial_mu_wrt_y(y,z); }
	fp partial_rho_wrt_z(fp x, fp y, fp z) const
		{ return local_coords::partial_mu_wrt_z(y,z); }
	fp partial_sigma_wrt_x(fp x, fp y, fp z) const
		{ return local_coords::partial_phi_wrt_x(x,y); }
	fp partial_sigma_wrt_y(fp x, fp y, fp z) const
		{ return local_coords::partial_phi_wrt_y(x,y); }
	fp partial_sigma_wrt_z(fp x, fp y, fp z) const { return 0.0; }

	// partial^2 (rho,sigma) / partial (xx,xy,xz,yy,yz)
	fp partial2_rho_wrt_xx(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_rho_wrt_xy(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_rho_wrt_xz(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_rho_wrt_yy(fp x, fp y, fp z) const
		{ return local_coords::partial2_mu_wrt_yy(y,z); }
	fp partial2_rho_wrt_yz(fp x, fp y, fp z) const
		{ return local_coords::partial2_mu_wrt_yz(y,z); }
	fp partial2_rho_wrt_zz(fp x, fp y, fp z) const
		{ return local_coords::partial2_mu_wrt_zz(y,z); }
	fp partial2_sigma_wrt_xx(fp x, fp y, fp z) const
		{ return local_coords::partial2_phi_wrt_xx(x,y); }
	fp partial2_sigma_wrt_xy(fp x, fp y, fp z) const
		{ return local_coords::partial2_phi_wrt_xy(x,y); }
	fp partial2_sigma_wrt_xz(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_sigma_wrt_yy(fp x, fp y, fp z) const
		{ return local_coords::partial2_phi_wrt_yy(x,y); }
	fp partial2_sigma_wrt_yz(fp x, fp y, fp z) const { return 0.0; }
	fp partial2_sigma_wrt_zz(fp x, fp y, fp z) const { return 0.0; }

	// plotting coordinates (px,py)
	// ... character string describing how (dpx,dpy) are
	//     defined in terms of (mu,nu,phi), eg "90 - drho = 90 - dphi"
	//     (used for labelling output files)
	const char* name_of_dpx() const
		{
		return is_plus() ? "90 - dsigma = 90 - dphi"
				 : "90 + dsigma = 90 + dphi";
		}
	const char* name_of_dpy() const { return "drho = dmu"; }
	// ... (rho,simga) --> (px,py)
	fp dpx_of_rho_sigma(fp rho, fp sigma) const
		{
		const fp dsigma = jtutil::degrees_of_radians(sigma);
		return is_plus() ? 90.0 - dsigma : 90.0 + dsigma;
		}
	fp dpy_of_rho_sigma(fp rho, fp sigma) const
		{ return jtutil::degrees_of_radians(rho); }

	// compute the arc length of a surface in the specified plane
	// (must be one of "xy" or "yz") over the patch's nominal bounds
	// ... error_exit() if  plane  is invalid
	fp plane_arc_length(const char plane[],
			    int ghosted_radius_gfn,
			    int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
					  int g_yy_gfn, int g_yz_gfn,
							int g_zz_gfn,
			    enum integration_method method)
		const;

	// constructor, destructor
	y_patch(patch_system &my_patch_system_in, int patch_number_in,
		const char* name_in, bool is_plus_in,
		const grid_arrays::grid_array_pars& grid_array_pars_in,
		const grid::grid_pars& grid_pars_in);
	~y_patch() { }

private:
        // we forbid copying and passing by value
        // by declaring the copy constructor and assignment operator
        // private, but never defining them
	y_patch(const y_patch& rhs);
	y_patch& operator=(const y_patch& rhs);
	};

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif
