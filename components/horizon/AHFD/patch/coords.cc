// coords.cc - coordinate conversions etc
// $Header$

//
// local_coords::fuzzy_EQ_{ang,dang}
// local_coords::modulo_reduce_{ang,dang}
//
// (r,(mu,nu,phi)) <--> (x,y,z)
// ((mu,nu,phi)) --> the 3rd
//
// partial_{x,y,z}_wrt_{mu,nu,phi}
// partial_{mu,nu,phi}_wrt_{x,y,z}
// partial2_{mu,nu,phi}_wrt_{xx,xy,xz,yy,yz,zz}
//
// usual polar spherical (r,theta,phi) <--> (x,y,z)
// ((mu,nu,phi)) <--> usual polar spherical (theta,phi)
// ((mu,nu,phi)) --> direction cosines (xcos,ycos,zcos)
//
// local_coords::mu_nu_phi::name_of_coords_set
//

#include <math.h>
#include <float.h>
#include <assert.h>
#include <limits.h>

#include "../cctk.h"

#include "../jtutil/util.hh"
using jtutil::arctan_xy;
using jtutil::signum;
using jtutil::pow2;
using jtutil::hypot3;

#include "coords.hh"

// all the code in this file is inside this namespace
namespace AHFD
	  {

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// these functions test if two angles are fuzzily equal mod 2*pi radians
// (360 degrees)
//

namespace local_coords
	  {

bool fuzzy_EQ_ang(fp ang1, fp ang2)
{
return jtutil::fuzzy<fp>::is_integer( (ang2 - ang1)/(2.0*PI) );
}

bool fuzzy_EQ_dang(fp dang1, fp dang2)
{
return jtutil::fuzzy<fp>::is_integer( (dang2 - dang1)/360.0 );
}

	  }

//******************************************************************************

//
// modulo-reduce  {ang,dang}  to be (fuzzily) within the range
// [min,max]_{ang,dang}, or error_exit() if no such value exists
//

namespace local_coords
	  {

fp modulo_reduce_ang(fp ang, fp min_ang, fp max_ang)
{
return jtutil::modulo_reduce(ang, 2.0*PI, min_ang, max_ang);
}

fp modulo_reduce_dang(fp dang, fp min_dang, fp max_dang)
{
return jtutil::modulo_reduce(dang, 360.0, min_dang, max_dang);
}

	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// these functions convert ((mu,nu,phi)) <--> (x,y,z)
//

//**************************************

// valid in +/- z patch
// not valid in xy plane (z == 0, i.e. (cos(mu) == 0 || cos(nu) == 0))
// unless r == 0 (x == y == z == 0)
// c.f. page 9.2 of my mpe notes
namespace local_coords
	  {
void xyz_of_r_mu_nu(fp r, fp mu, fp nu,   fp& x, fp& y, fp& z)
{
const fp sign_y        = signum(sin(mu));
const fp sign_z_via_mu = signum(cos(mu));
assert( jtutil::fuzzy<fp>::NE(cos(mu), 0.0) );
const fp y_over_z = tan(mu);

const fp sign_x        = signum(sin(nu));
const fp sign_z_via_nu = signum(cos(nu));
assert( jtutil::fuzzy<fp>::NE(cos(nu), 0.0) );
const fp x_over_z = tan(nu);

// failure of next assert() ==> inconsistent input (mu,nu)
assert(sign_z_via_mu == sign_z_via_nu);
const fp sign_z = sign_z_via_mu;

const fp temp = 1.0 / sqrt(1.0 + pow2(y_over_z) + pow2(x_over_z));

z = sign_z * r * temp;
x = x_over_z * z;
y = y_over_z * z;

// if (jtutil::fuzzy<fp>::NE(r, 0.0))
//    then {
// 	if (jtutil::fuzzy<fp>::NE(x, 0.0))
//  	   then assert( signum(x) == sign_x );
// 	if (jtutil::fuzzy<fp>::NE(y, 0.0))
// 	   then assert( signum(y) == sign_y );
// 	}
}
	  }

//**************************************

// valid in +/- y patch
// not valid in xz plane (y == 0, i.e. (sin(mu) == 0 || sin(phi) == 0)
// unless r == 0 (x == y == z == 0)
// c.f. page 9.3 of my mpe notes
namespace local_coords
	  {
void xyz_of_r_mu_phi(fp r, fp mu, fp phi,   fp& x, fp& y, fp& z)
{
const fp  mu_bar = 0.5*PI - mu ;
const fp phi_bar = 0.5*PI - phi;

const fp sign_z            = signum(sin(mu_bar));
const fp sign_y_via_mu_bar = signum(cos(mu_bar));
assert( jtutil::fuzzy<fp>::NE(cos(mu_bar), 0.0) );
const fp z_over_y = tan(mu_bar);

const fp sign_x             = signum(sin(phi_bar));
const fp sign_y_via_phi_bar = signum(cos(phi_bar));
assert( jtutil::fuzzy<fp>::NE(cos(phi_bar), 0.0) );
const fp x_over_y = tan(phi_bar);

// failure of next assert() ==> inconsistent input (mu,phi)
assert(sign_y_via_mu_bar == sign_y_via_phi_bar);
const fp sign_y = sign_y_via_mu_bar;

const fp temp = 1.0 / sqrt(1.0 + pow2(z_over_y) + pow2(x_over_y));

y = sign_y * r * temp;
z = z_over_y * y;
x = x_over_y * y;

// if (jtutil::fuzzy<fp>::NE(r, 0.0))
//    then {
// 	if (jtutil::fuzzy<fp>::NE(z, 0.0))
// 	   then assert( signum(z) == sign_z );
// 	if (jtutil::fuzzy<fp>::NE(x, 0.0))
// 	   then assert( signum(x) == sign_x );
// 	}
}
	  }

//**************************************

// valid in +/- x patch
// not valid in yz plane (x == 0, i.e. (sin(nu) == 0 || cos(phi) == 0))
// unless r == 0 (x == y == z == 0)
// c.f. page 9.4 of my mpe notes
namespace local_coords
	  {
void xyz_of_r_nu_phi(fp r, fp nu, fp phi,   fp& x, fp& y, fp& z)
{
const fp nu_bar = 0.5*PI - nu ;

const fp sign_z            = signum(sin(nu_bar));
const fp sign_x_via_nu_bar = signum(cos(nu_bar));
assert( jtutil::fuzzy<fp>::NE(cos(nu_bar), 0.0) );
const fp z_over_x = tan(nu_bar);

const fp sign_y         = signum(sin(phi));
const fp sign_x_via_phi = signum(cos(phi));
assert( jtutil::fuzzy<fp>::NE(cos(phi), 0.0) );
const fp y_over_x = tan(phi   );

// failure of next assert() ==> inconsistent input (nu,phi)
assert(sign_x_via_nu_bar == sign_x_via_phi);
const fp sign_x = sign_x_via_nu_bar;

const fp temp = 1.0 / sqrt(1.0 + pow2(z_over_x) + pow2(y_over_x));

x = sign_x * r * temp;
z = z_over_x * x;
y = y_over_x * x;

// if (jtutil::fuzzy<fp>::NE(r, 0.0))
//    then {
// 	if (jtutil::fuzzy<fp>::NE(z, 0.0))
// 	   then assert( signum(z) == sign_z );
// 	if (jtutil::fuzzy<fp>::NE(y, 0.0))
// 	   then assert( signum(y) == sign_y );
// 	}
}
	  }

//******************************************************************************

//
// these functions give any one of {mu,nu,phi} in terms of the other two
//

//**************************************

// ill-conditioned near z axis
// not valid in xy plane (cos(mu) == 0 || cos(nu) == 0)
namespace local_coords
	  {
fp phi_of_mu_nu(fp mu, fp nu)
{
fp x, y, z;
xyz_of_r_mu_nu(1.0,mu,nu, x,y,z);
return phi_of_xy(x, y);
}
	  }

//**************************************

// ill-conditioned near y axis
// not valid in xz plane (sin(mu) == 0 || sin(phi == 0)
namespace local_coords
	  {
fp nu_of_mu_phi(fp mu, fp phi)
{
fp x, y, z;
xyz_of_r_mu_phi(1.0,mu,phi, x,y,z);
return nu_of_xz(x, z);
}
	  }

//**************************************

// ill-conditioned near x axis
// not valid in yz plane (sin(nu) == 0 || cos(phi) == 0)
namespace local_coords
	  {
fp mu_of_nu_phi(fp nu, fp phi)
{
fp x, y, z;
xyz_of_r_nu_phi(1.0,nu,phi, x,y,z);
return mu_of_yz(y, z);
}
	  }


//******************************************************************************

namespace local_coords
	  {
fp r_of_xyz(fp x, fp y, fp z) { return hypot3(x, y, z); }
fp mu_of_yz(fp y, fp z)  { return arctan_xy(z,y); }
fp nu_of_xz(fp x, fp z)  { return arctan_xy(z,x); }
fp phi_of_xy(fp x, fp y) { return arctan_xy(x,y); }
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// these functions compute the partial derivatives
//	partial {x,y,z} / partial {mu,nu,phi}
//
// Bugs: they're slow :(
//

namespace local_coords
	  {
void partial_xyz_wrt_r_mu_nu
	(fp r, fp mu, fp nu,
	 fp& partial_x_wrt_r, fp& partial_x_wrt_mu, fp& partial_x_wrt_nu,
	 fp& partial_y_wrt_r, fp& partial_y_wrt_mu, fp& partial_y_wrt_nu,
	 fp& partial_z_wrt_r, fp& partial_z_wrt_mu, fp& partial_z_wrt_nu)
{
const fp tan_mu = tan(mu);
const fp tan_nu = tan(nu);
const fp tan2_mu = pow2(tan_mu);
const fp tan2_nu = pow2(tan_nu);

fp x, y, z;
xyz_of_r_mu_nu(r,mu,nu, x,y,z);

// Comment this out, accept a nan instead
// assert( jtutil::fuzzy<fp>::NE(r, 0.0) );
const fp rinv = 1.0/r;
partial_x_wrt_r = x*rinv;
partial_y_wrt_r = y*rinv;
partial_z_wrt_r = z*rinv;

const fp t = 1 + tan2_mu + tan2_nu;	// = $r^2/z^2$
const fp partial_t_wrt_mu = 2.0 * tan_mu * (1.0+tan2_mu);
const fp partial_t_wrt_nu = 2.0 * tan_nu * (1.0+tan2_nu);

const fp r2_over_zt2 = (r*r) / (z*t*t);
partial_z_wrt_mu = -0.5 * r2_over_zt2 * partial_t_wrt_mu;
partial_z_wrt_nu = -0.5 * r2_over_zt2 * partial_t_wrt_nu;

partial_x_wrt_mu = tan_nu*partial_z_wrt_mu;
partial_x_wrt_nu = tan_nu*partial_z_wrt_nu + z*(1.0+tan2_nu);
partial_y_wrt_mu = tan_mu*partial_z_wrt_mu + z*(1.0+tan2_mu);
partial_y_wrt_nu = tan_mu*partial_z_wrt_nu;
}
	  }

//**************************************

namespace local_coords
	  {
void partial_xyz_wrt_r_mu_phi
	(fp r, fp mu, fp phi,
	 fp& partial_x_wrt_r, fp& partial_x_wrt_mu, fp& partial_x_wrt_phi,
	 fp& partial_y_wrt_r, fp& partial_y_wrt_mu, fp& partial_y_wrt_phi,
	 fp& partial_z_wrt_r, fp& partial_z_wrt_mu, fp& partial_z_wrt_phi)
{
const fp  mu_bar = 0.5*PI - mu ;
const fp phi_bar = 0.5*PI - phi;

const fp tan_mu_bar  = tan(mu_bar);
const fp tan_phi_bar = tan(phi_bar);
const fp tan2_mu_bar  = pow2(tan_mu_bar);
const fp tan2_phi_bar = pow2(tan_phi_bar);

fp x, y, z;
xyz_of_r_mu_phi(r,mu,phi, x,y,z);

// Comment this out, accept a nan instead
// assert( jtutil::fuzzy<fp>::NE(r, 0.0) );
const fp rinv = 1.0/r;
partial_x_wrt_r = x*rinv;
partial_y_wrt_r = y*rinv;
partial_z_wrt_r = z*rinv;

const fp t = 1 + tan2_mu_bar + tan2_phi_bar;	// = $r^2/y^2$
const fp partial_t_wrt_mu_bar  = 2.0 * tan_mu_bar  * (1.0+tan2_mu_bar );
const fp partial_t_wrt_phi_bar = 2.0 * tan_phi_bar * (1.0+tan2_phi_bar);

const fp r2_over_yt2 = (r*r) / (y*t*t);
partial_y_wrt_mu  = 0.5 * r2_over_yt2 * partial_t_wrt_mu_bar;
partial_y_wrt_phi = 0.5 * r2_over_yt2 * partial_t_wrt_phi_bar;

partial_x_wrt_mu  = tan_phi_bar*partial_y_wrt_mu;
partial_x_wrt_phi = tan_phi_bar*partial_y_wrt_phi - y*(1.0+tan2_phi_bar);
partial_z_wrt_mu  = tan_mu_bar *partial_y_wrt_mu  - y*(1.0+tan2_mu_bar );
partial_z_wrt_phi = tan_mu_bar *partial_y_wrt_phi;
}
	  }

//**************************************

namespace local_coords
	  {
void partial_xyz_wrt_r_nu_phi
	(fp r, fp nu, fp phi,
	 fp& partial_x_wrt_r, fp& partial_x_wrt_nu, fp& partial_x_wrt_phi,
	 fp& partial_y_wrt_r, fp& partial_y_wrt_nu, fp& partial_y_wrt_phi,
	 fp& partial_z_wrt_r, fp& partial_z_wrt_nu, fp& partial_z_wrt_phi)
{
const fp  nu_bar = 0.5*PI - nu ;

const fp tan_nu_bar = tan(nu_bar);
const fp tan_phi    = tan(phi);
const fp tan2_nu_bar = pow2(tan_nu_bar);
const fp tan2_phi    = pow2(tan_phi);

fp x, y, z;
xyz_of_r_nu_phi(r,nu,phi, x,y,z);

// Comment this out, accept a nan instead
// assert( jtutil::fuzzy<fp>::NE(r, 0.0) );
const fp rinv = 1.0/r;
partial_x_wrt_r = x*rinv;
partial_y_wrt_r = y*rinv;
partial_z_wrt_r = z*rinv;

const fp t = 1 + tan2_nu_bar + tan2_phi;	// = $r^2/x^2$
const fp partial_t_wrt_nu_bar  = 2.0 * tan_nu_bar * (1.0+tan2_nu_bar);
const fp partial_t_wrt_phi     = 2.0 * tan_phi    * (1.0+tan2_phi   );

const fp r2_over_xt2 = (r*r) / (x*t*t);
partial_x_wrt_nu  =  0.5 * r2_over_xt2 * partial_t_wrt_nu_bar;
partial_x_wrt_phi = -0.5 * r2_over_xt2 * partial_t_wrt_phi;

partial_y_wrt_nu  = tan_phi   *partial_x_wrt_nu;
partial_y_wrt_phi = tan_phi   *partial_x_wrt_phi + x*(1.0+tan2_phi   );
partial_z_wrt_nu  = tan_nu_bar*partial_x_wrt_nu  - x*(1.0+tan2_nu_bar);
partial_z_wrt_phi = tan_nu_bar*partial_x_wrt_phi;
}
	  }

//******************************************************************************

//
// these functions compute the partial derivatives
//	partial {mu,nu,phi} / partial {x,y,z}
// as computed by the maple file "coord_derivs.{maple,out}" in this directory
//
namespace local_coords
	  {
fp partial_mu_wrt_y(fp y, fp z) { return  z / (y*y + z*z); }
fp partial_mu_wrt_z(fp y, fp z) { return -y / (y*y + z*z); }

fp partial_nu_wrt_x(fp x, fp z) { return  z / (x*x + z*z); }
fp partial_nu_wrt_z(fp x, fp z) { return -x / (x*x + z*z); }

fp partial_phi_wrt_x(fp x, fp y) { return -y / (x*x + y*y); }
fp partial_phi_wrt_y(fp x, fp y) { return  x / (x*x + y*y); }
	  }

//******************************************************************************

//
// these functions compute the 2nd partial derivatives
//	partial {mu,nu,phi} / partial {xx,xy,xz,yy,yz,zz}
// as computed by the maple file "coord_derivs.{maple,out}" in this directory
//
namespace local_coords
	  {
fp partial2_mu_wrt_yy(fp y, fp z) {return   -2.0*y*z  / pow2(y*y + z*z);}
fp partial2_mu_wrt_yz(fp y, fp z) {return (y*y - z*z) / pow2(y*y + z*z);}
fp partial2_mu_wrt_zz(fp y, fp z) {return    2.0*y*z  / pow2(y*y + z*z);}

fp partial2_nu_wrt_xx(fp x, fp z) {return   -2.0*x*z  / pow2(x*x + z*z);}
fp partial2_nu_wrt_xz(fp x, fp z) {return (x*x - z*z) / pow2(x*x + z*z);}
fp partial2_nu_wrt_zz(fp x, fp z) {return    2.0*x*z  / pow2(x*x + z*z);}

fp partial2_phi_wrt_xx(fp x, fp y) {return    2.0*x*y  / pow2(x*x + y*y);}
fp partial2_phi_wrt_xy(fp x, fp y) {return (y*y - x*x) / pow2(x*x + y*y);}
fp partial2_phi_wrt_yy(fp x, fp y) {return   -2.0*x*y  / pow2(x*x + y*y);}
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// these functions convert the usual polar spherical (theta,phi) <--> (x,y,z)
//

//**************************************

namespace local_coords
	  {
void xyz_of_r_theta_phi(fp r, fp theta, fp phi,   fp& x, fp& y, fp& z)
{
z = r * cos(theta);
x = r * sin(theta) * cos(phi);
y = r * sin(theta) * sin(phi);
}
	  }

//**************************************

namespace local_coords
	  {
void r_theta_phi_of_xyz(fp x, fp y, fp z,   fp& r, fp& theta, fp& phi)
{
r = r_of_xyz(x, y, z);
theta = theta_of_xyz(x, y, z);
phi = phi_of_xy(x, y);
}
	  }

//**************************************

namespace local_coords
	  {
fp theta_of_xyz(fp x, fp y, fp z)
{
return arctan_xy(z , hypot(x,y));
}
	  }

//******************************************************************************

//
// these functions convert ((mu,nu,phi)) <--> usual polar spherical (theta,phi)
// ... note phi is the same coordinate in both systems
//

namespace local_coords
	  {
void theta_phi_of_mu_nu(fp mu, fp nu,   fp& ps_theta, fp& ps_phi)
{
fp x, y, z;
xyz_of_r_mu_nu(1.0,mu,nu, x,y,z);

ps_theta = theta_of_xyz(x, y, z);
ps_phi = phi_of_xy(x, y);
}
	  }

//**************************************

// Bugs: computes ps_phi via trig, even though it's trivially == phi
namespace local_coords
	  {
void theta_phi_of_mu_phi(fp mu, fp phi,   fp& ps_theta, fp& ps_phi)
{
fp x, y, z;
xyz_of_r_mu_phi(1.0,mu,phi, x,y,z);

ps_theta = theta_of_xyz(x, y, z);
ps_phi = phi_of_xy(x, y);
assert( fuzzy_EQ_ang(ps_phi, phi) );
}
	  }

//**************************************

// Bugs: computes ps_phi via trig, even though it's trivially == phi
namespace local_coords
	  {
void theta_phi_of_nu_phi(fp nu, fp phi,   fp& ps_theta, fp& ps_phi)
{
fp x, y, z;
xyz_of_r_nu_phi(1.0,nu,phi, x,y,z);

ps_theta = theta_of_xyz(x, y, z);
ps_phi = phi_of_xy(x, y);
assert( fuzzy_EQ_ang(ps_phi, phi) );
}
	  }

//******************************************************************************

namespace local_coords
	  {
void mu_nu_of_theta_phi(fp ps_theta, fp ps_phi,   fp& mu, fp& nu)
{
fp x, y, z;
xyz_of_r_theta_phi(1.0,ps_theta,ps_phi, x,y,z);

mu = mu_of_yz(y, z);
nu = nu_of_xz(x, z);
}
	  }

//**************************************

// Bugs: computes phi via trig, even though it's trivially == ps_phi
namespace local_coords
	  {
void mu_phi_of_theta_phi(fp ps_theta, fp ps_phi,   fp& mu, fp& phi)
{
fp x, y, z;
xyz_of_r_theta_phi(1.0,ps_theta,ps_phi, x,y,z);

mu  =  mu_of_yz(y, z);
phi = phi_of_xy(x, y);
assert( fuzzy_EQ_ang(phi, ps_phi) );
}
	  }

//**************************************

// Bugs: computes phi via trig, even though it's trivially == ps_phi
namespace local_coords
	  {
void nu_phi_of_theta_phi(fp ps_theta, fp ps_phi,   fp& nu, fp& phi)
{
fp x, y, z;
xyz_of_r_theta_phi(1.0,ps_theta,ps_phi, x,y,z);

nu  =  nu_of_xz(x, z);
phi = phi_of_xy(x, y);
assert( fuzzy_EQ_ang(phi, ps_phi) );
}
	  }

//******************************************************************************

//
// these functions convert ((mu,nu,phi)) to the direction cosines
// (xcos,ycos,zcos)
//

namespace local_coords
	  {
void xyzcos_of_mu_nu (fp mu, fp nu , fp& xcos, fp& ycos, fp& zcos)
{
xyz_of_r_mu_nu(1.0,mu,nu, xcos,ycos,zcos);
}
	  }

namespace local_coords
	  {
void xyzcos_of_mu_phi(fp mu, fp phi, fp& xcos, fp& ycos, fp& zcos)
{
xyz_of_r_mu_phi(1.0,mu,phi, xcos,ycos,zcos);
}
	  }

namespace local_coords
	  {
void xyzcos_of_nu_phi(fp nu, fp phi, fp& xcos, fp& ycos, fp& zcos)
{
xyz_of_r_nu_phi(1.0,nu,phi, xcos,ycos,zcos);
}
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function computes a human-readable name from a (mu,nu,phi)
// coordinates set.
//
const char* local_coords::name_of_coords_set(coords_set S)
{
//
// we have to use an if-else chain because the  local_coords::set_*
// constants aren't compile-time constants and hence aren't eligible
// to be switch case labels
//
if	(S == coords_set_empty)
   then return "{}";
else if (S == coords_set_mu)
   then return "mu";
else if (S == coords_set_nu)
   then return "nu";
else if (S == coords_set_phi)
   then return "phi";
else if (S == (coords_set_mu|coords_set_nu))
   then return "{mu,nu}";
else if (S == (coords_set_mu|coords_set_phi))
   then return "{mu,phi}";
else if (S == (coords_set_nu|coords_set_phi))
   then return "{nu,phi}";
else if (S == (coords_set_mu|coords_set_nu|coords_set_phi))
   then return "{mu,nu,phi}";
else	error_exit(PANIC_EXIT,
"***** local_coords::mu_nu_phi::name_of_coords_set:\n"
"        S=0x%x isn't a valid  coords_set  bit vector!\n"
,
		   int(S));					/*NOTREACHED*/
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

	  }	// namespace AHFinderDirect
