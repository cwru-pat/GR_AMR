#ifndef AHFD_PATCH_COORDS_H
#define AHFD_PATCH_COORDS_H
// coords.hh -- coordinate systems and conversions
// $Header$
//
// ***** coordinate systems *****
// ***** conversions between different local coordinate systems *****
// ***** bit masks for coordinates ****
// global_coords - conversions between global and local coordinates
//

//
// prerequisites:
//    <math.h>
//    "cctk.h" or "fake_cctk.h"
//    "stdc.h"
//    "config.hh"
//    "../jtutil/util.hh"
//

// everything in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

//
// ***** coordinate systems *****
//

//
// We use the following terminology for coordinates:
//	{a,b,c} == any one of (a,b,c) == a | b | c
//	((a,b,c)) == any two of (a,b,c) == (a,b) | (a,c) | (b,c)
//	(r,(a,b,c)) == r and any two of (a,b,c) == (r,a,b) | (r,a,c) | (r,b,c)
//

//
// We use the following coordinates and coordinate systems:
//	global_(x,y,z)	global Cartesian x,y,z coordinates
//			(these are the coordinates Cactus gives us)
//	(x,y,z)		local Cartesian x,y,z coordinates relative to a
//			local origin point (stored in this object);
//			This origin should be chosen to be somewhere
//			at least vaguely in the "middle" of our
//			apparent horizon
//	(r,theta,phi)	polar spherical coordinates about the local
//			origin point
//	(r,(mu,nu,phi))	patch coordinates about the local origin point
//		mu	arctan(y/z) = rotation about local x axis
//		nu	arctan(x/z) = rotation about local y axis
//		phi	arctan(y/x) = rotation about local z axis
//				    = usual polar spherical phi
//	(rho,sigma)	generic patch angular coordinates, these will
//			be whichever two of (mu,nu,phi) are nonsingular
//			throughout the patch
//	tau		the third member of (mu,nu,phi) (the one which
//			is neither rho nor sigma)
//	(perp,par)	(rho,sigma) coordinates (perpendicular,parallel)
//			to a patch boundary; these are defined by class
//			patch_edge and used by a lot of the boundary code
//

//
// We also use direction cosines (xcos,ycos,zcos):
// xcos = x/r = cos(angle between vector from origin to (x,y,z), and x axis)
// ycos = y/r = cos(angle between vector from origin to (x,y,z), and y axis)
// zcos = z/r = cos(angle between vector from origin to (x,y,z), and z axis)
//

//
// We also use prefixes on coordinates, eg.
//	irho		integer grid coordinate
//	 rho		floating point angular coordinate, measured in radians
//	drho		floating point angular coordinate, measured in degrees
// I/O uses the degree coordinates in preference to radian ones.
//

//
// We also define plotting coordinates (px,py) (in both fp and degrees-fp
// variants), based on the conventional "looking down the z axis" picture
// of a patch system.  Unfortunately there's no good place to stick a -z
// patch in this scheme, so we arbitrarily put it just beyond the +x patch:
//
//	           py axis
//	              |
//	          +-------+
//	          |       |
//	          |  +y   |
//	          |       |
//	  +-------+-------+-------+-------+
//	  |       |       |       |       |
//	--|  -x   |  +z   |  +x   |  -z   |-- px axis
//	  |       |       |       |       |
//	  +-------+-------+-------+-------+
//	          |       |
//	          |  -y   |
//	          |       |
//	          +-------+
//	              |
//
// The main (in fact only) use of these coordinates is in their degrees
// forms (dpx,dpy), for I/O and constructor arguments.
//

//*****************************************************************************

//
// ***** stuff for treating angles modulo 2*pi radians (360 degrees) *****
//

namespace local_coords
	  {

// compare if two angles are fuzzily equal mod 2*pi radians (360 degrees)
bool fuzzy_EQ_ang(fp ang1, fp ang2);	// radians
bool fuzzy_EQ_dang(fp dang1, fp dang2);	// degrees

// modulo-reduce  {ang,dang}  to be (fuzzily) within the range
// [min,max]_{ang,dang}, or error_exit() if no such value exists
fp modulo_reduce_ang(fp ang, fp min_ang, fp max_ang);
fp modulo_reduce_dang(fp dang, fp min_dang, fp max_dang);

	  }	// close namespace local_coords::

//*****************************************************************************

//
// ***** conversions between different local coordinate systems *****
//
// Bugs:
//	These functions aren't optimally efficient: they do lots
//	of could-be-optimized away arithmetic and trig calls.  But
//	In practice this doesn't matter, since we don't call these
//	functions inside inner loops (instead, we cache their values
//	at grid points).
//

namespace local_coords
	  {
// (r,(mu,nu,phi)) <--> (x,y,z)
void xyz_of_r_mu_nu (fp r, fp mu, fp nu ,   fp& x, fp& y, fp& z);
void xyz_of_r_mu_phi(fp r, fp mu, fp phi,   fp& x, fp& y, fp& z);
void xyz_of_r_nu_phi(fp r, fp nu, fp phi,   fp& x, fp& y, fp& z);
fp r_of_xyz(fp x, fp y, fp z);
fp mu_of_yz(fp y, fp z);
fp nu_of_xz(fp x, fp z);
fp phi_of_xy(fp x, fp y);

// ((mu,nu,phi)) --> the 3rd
fp phi_of_mu_nu(fp mu, fp nu );
fp nu_of_mu_phi(fp mu, fp phi);
fp mu_of_nu_phi(fp nu, fp phi);

// partial {x,y,z} / partial {mu,nu,phi}
void partial_xyz_wrt_r_mu_nu
	(fp r, fp mu, fp nu,
	 fp& partial_x_wrt_r, fp& partial_x_wrt_mu, fp& partial_x_wrt_nu,
	 fp& partial_y_wrt_r, fp& partial_y_wrt_mu, fp& partial_y_wrt_nu,
	 fp& partial_z_wrt_r, fp& partial_z_wrt_mu, fp& partial_z_wrt_nu);
void partial_xyz_wrt_r_mu_phi
	(fp r, fp mu, fp phi,
	 fp& partial_x_wrt_r, fp& partial_x_wrt_mu, fp& partial_x_wrt_phi,
	 fp& partial_y_wrt_r, fp& partial_y_wrt_mu, fp& partial_y_wrt_phi,
	 fp& partial_z_wrt_r, fp& partial_z_wrt_mu, fp& partial_z_wrt_phi);
void partial_xyz_wrt_r_nu_phi
	(fp r, fp nu, fp phi,
	 fp& partial_x_wrt_r, fp& partial_x_wrt_nu, fp& partial_x_wrt_phi,
	 fp& partial_y_wrt_r, fp& partial_y_wrt_nu, fp& partial_y_wrt_phi,
	 fp& partial_z_wrt_r, fp& partial_z_wrt_nu, fp& partial_z_wrt_phi);

// partial {mu,nu,phi} / partial {x,y,z}
fp partial_mu_wrt_y(fp y, fp z);
fp partial_mu_wrt_z(fp y, fp z);
fp partial_nu_wrt_x(fp x, fp z);
fp partial_nu_wrt_z(fp x, fp z);
fp partial_phi_wrt_x(fp x, fp y);
fp partial_phi_wrt_y(fp x, fp y);

// partial^2 {mu,nu,phi} / partial {x,y,z}{x,y,z}
fp partial2_mu_wrt_yy(fp y, fp z);
fp partial2_mu_wrt_yz(fp y, fp z);
fp partial2_mu_wrt_zz(fp y, fp z);
fp partial2_nu_wrt_xx(fp x, fp z);
fp partial2_nu_wrt_xz(fp x, fp z);
fp partial2_nu_wrt_zz(fp x, fp z);
fp partial2_phi_wrt_xx(fp x, fp y);
fp partial2_phi_wrt_xy(fp x, fp y);
fp partial2_phi_wrt_yy(fp x, fp y);

// usual polar spherical (r,theta,phi) <--> (x,y,z)
void xyz_of_r_theta_phi(fp r, fp theta, fp phi,   fp& x, fp& y, fp& z);
void r_theta_phi_of_xyz(fp x, fp y, fp z,   fp& r, fp& theta, fp& phi);
// ... already have r_of_xyz()
// ... already have phi_of_xy()
fp theta_of_xyz(fp x, fp y, fp z);

// ((mu,nu,phi)) <--> usual polar spherical (theta,phi)
// ... note phi is the same coordinate in both systems
void theta_phi_of_mu_nu (fp mu, fp nu ,   fp& ps_theta, fp& ps_phi);
void theta_phi_of_mu_phi(fp mu, fp phi,   fp& ps_theta, fp& ps_phi);
void theta_phi_of_nu_phi(fp nu, fp phi,   fp& ps_theta, fp& ps_phi);
void  mu_nu_of_theta_phi(fp ps_theta, fp ps_phi,   fp& mu, fp& nu );
void mu_phi_of_theta_phi(fp ps_theta, fp ps_phi,   fp& mu, fp& phi);
void nu_phi_of_theta_phi(fp ps_theta, fp ps_phi,   fp& nu, fp& phi);

// ((mu,nu,phi)) --> direction cosines (xcos,ycos,zcos)
void xyzcos_of_mu_nu (fp mu, fp nu , fp& xcos, fp& ycos, fp& zcos);
void xyzcos_of_mu_phi(fp mu, fp phi, fp& xcos, fp& ycos, fp& zcos);
void xyzcos_of_nu_phi(fp nu, fp phi, fp& xcos, fp& ycos, fp& zcos);
	  }	// close namespace local_coords::

//*****************************************************************************

//
// ***** bit masks for coordinates ****
//

//
// We need to manipulate coordinates to do calculations like "which
// coordinate do these two patches have in common".  We do these by
// Boolean operations on integers using the following bit masks:
//

namespace local_coords
	  {

typedef int coords_set;

enum {
     coords_set_mu  = 0x1,
     coords_set_nu  = 0x2,
     coords_set_phi = 0x4,

     coords_set_empty = 0x0,
     coords_set_all = coords_set_mu | coords_set_nu | coords_set_phi // no comma
     };

// human-readable coordinate names for debugging etc
const char* name_of_coords_set(coords_set S);

// set complement of coordinates
inline
  coords_set coords_set_not(coords_set S)
	{ return coords_set_all & ~S; }

	  }	// close namespace local_coords::

//******************************************************************************

//
// This class stores the origin point of our local coordinates, and
// provides conversions between local and global coordinates.
//
class	global_coords
	{
public:
	#ifdef NOT_USED
	// global (x,y,z) --> local (x,y,z)
	fp local_x_of_global_x(fp global_x) const
		{ return global_x - origin_x_; }
	fp local_y_of_global_y(fp global_y) const
		{ return global_y - origin_y_; }
	fp local_z_of_global_z(fp global_z) const
		{ return global_z - origin_z_; }
	#endif	/* NOT_USED */

	#ifdef NOT_USED
	// local (x,y,z) --> global (x,y,z)
	fp global_x_of_local_x(fp local_x) const
		{ return origin_x_ + local_x; }
	fp global_y_of_local_y(fp local_y) const
		{ return origin_y_ + local_y; }
	fp global_z_of_local_z(fp local_z) const
		{ return origin_z_ + local_z; }
	#endif	/* NOT_USED */

	// get global (x,y,z) coordinates of local origin point
	fp origin_x() const { return origin_x_; }
	fp origin_y() const { return origin_y_; }
	fp origin_z() const { return origin_z_; }

	// set global (x,y,z) coordinates of local origin point
	void origin_x(const fp x) { origin_x_=x; }
	void origin_y(const fp y) { origin_y_=y; }
	void origin_z(const fp z) { origin_z_=z; }

	// radius of given (x,y,z) with respect to local origin point
	#ifdef NOT_USED
	fp radius_of_local_xyz(fp local_x, fp local_y, fp local_z) const
		{ return jtutil::hypot3(local_x, local_y, local_z); }
	fp radius_of_global_xyz(fp global_x, fp global_y, fp global_z)
		const
		{
		return radius_of_local_xyz(local_x_of_global_x(global_x),
					   local_y_of_global_y(global_y),
					   local_z_of_global_z(global_z));
		}
	#endif	/* NOT_USED */

	// constructor: specify global (x,y,z) coordinates of local origin point
	global_coords(fp origin_x_in, fp origin_y_in, fp origin_z_in)
		: origin_x_(origin_x_in),
		  origin_y_(origin_y_in),
		  origin_z_(origin_z_in)
		{ }
	// destructor: compiler-generated no-op is ok

private:
	// we forbid copying and passing by value
	// by declaring the copy constructor and assignment operator
	// private, but never defining them
	global_coords(const global_coords& rhs);
	global_coords& operator=(const global_coords& rhs);

private:
	// global (x,y,z) coordinates of local origin point
	fp origin_x_, origin_y_, origin_z_;
	};

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif
