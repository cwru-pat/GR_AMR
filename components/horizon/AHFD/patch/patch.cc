// patch.cc -- describes a coordinate/grid patch
// $Header$

//
// patch::patch
// patch::~patch
// z_patch::z_patch
// x_patch::x_patch
// y_patch::y_patch
//
// patch::rho_sigma_metric
//
// patch::decode_integration_method
// patch::rho_arc_length
// patch::sigma_arc_length
// z_patch::plane_arc_length
// x_patch::plane_arc_length
// y_patch::plane_arc_length
// patch::integrate_gridfn
/// patch::integration_coeff
//
// patch::ghost_zone_on_edge
// patch::corner_ghost_zone_containing_point
// patch::ghost_zone_containing_noncorner_point
// patch::create_mirror_symmetry_ghost_zone
// patch::create_periodic_symmetry_ghost_zone
// patch::create_interpatch_ghost_zone
/// patch::set_ghost_zone
// patch::edge_adjacent_to_patch
// patch::assert_all_ghost_zones_fully_setup
//

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

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

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function constructs a  patch  object.
//
patch::patch(patch_system &my_patch_system_in, int patch_number_in,
	     const char name_in[], bool is_plus_in, char ctype_in,
	     local_coords::coords_set coords_set_rho_in,
	     local_coords::coords_set coords_set_sigma_in,
	     local_coords::coords_set coords_set_tau_in,
	     const grid_arrays::grid_array_pars& grid_array_pars_in,
	     const grid::grid_pars& grid_pars_in)

	: fd_grid(grid_array_pars_in, grid_pars_in),

	  my_patch_system_(my_patch_system_in),
	  patch_number_(patch_number_in),
	  name_(name_in),
	  is_plus_(is_plus_in), ctype_(ctype_in),

	  coords_set_rho_  (coords_set_rho_in  ),
	  coords_set_sigma_(coords_set_sigma_in),
	  coords_set_tau_  (coords_set_tau_in  ),

	  min_rho_patch_edge_(*new patch_edge(*this, side_is_min, side_is_rho)),
	  max_rho_patch_edge_(*new patch_edge(*this, side_is_max, side_is_rho)),
	  min_sigma_patch_edge_
	  	(*new patch_edge(*this, side_is_min, side_is_sigma)),
	  max_sigma_patch_edge_
		(*new patch_edge(*this, side_is_max, side_is_sigma)),

	  min_rho_ghost_zone_(NULL),
	  max_rho_ghost_zone_(NULL),
	  min_sigma_ghost_zone_(NULL),
	  max_sigma_ghost_zone_(NULL) // no comma

{ }

//******************************************************************************

//
// This function destroys a  patch  object.
//
patch::~patch()
{
// no need to check for null pointers, since  delete NULL  is a silent no-op

delete max_sigma_ghost_zone_;
delete min_sigma_ghost_zone_;
delete max_rho_ghost_zone_;
delete min_rho_ghost_zone_;

delete & max_sigma_patch_edge_;
delete & min_sigma_patch_edge_;
delete & max_rho_patch_edge_;
delete & min_rho_patch_edge_;
}

//******************************************************************************

//
// This function constructs a  z_patch  object.
//
z_patch::z_patch(patch_system &my_patch_system_in, int patch_number_in,
		 const char* name_in, bool is_plus_in,
		 const grid_arrays::grid_array_pars& grid_array_pars_in,
		 const grid::grid_pars& grid_pars_in)
	: patch(my_patch_system_in, patch_number_in,
		name_in, is_plus_in, 'z',
		local_coords::coords_set_mu, local_coords::coords_set_nu,
		local_coords::coords_set_phi,
		grid_array_pars_in, grid_pars_in)
{ }

//******************************************************************************

//
// This function constructs an  x_patch  object.
//
x_patch::x_patch(patch_system &my_patch_system_in, int patch_number_in,
		 const char* name_in, bool is_plus_in,
		 const grid_arrays::grid_array_pars& grid_array_pars_in,
		 const grid::grid_pars& grid_pars_in)
	: patch(my_patch_system_in, patch_number_in,
		name_in, is_plus_in, 'x',
		local_coords::coords_set_nu, local_coords::coords_set_phi,
		local_coords::coords_set_mu,
		grid_array_pars_in, grid_pars_in)
{ }

//******************************************************************************

//
// This function constructs a  y_patch  object.
//
y_patch::y_patch(patch_system &my_patch_system_in, int patch_number_in,
		 const char* name_in, bool is_plus_in,
		 const grid_arrays::grid_array_pars& grid_array_pars_in,
		 const grid::grid_pars& grid_pars_in)
	: patch(my_patch_system_in, patch_number_in,
		name_in, is_plus_in, 'y',
		local_coords::coords_set_mu, local_coords::coords_set_phi,
		local_coords::coords_set_nu,
		grid_array_pars_in, grid_pars_in)
{ }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function computes the (rho,sigma) induced 2-D metric from the
// 3-D (x,y,z) metric of the space containing the patch, as per p.33 of
// my apparent horizon finding notes.
//
// Arguments:
// (r,rho,sigma) = The coordinates where the Jacobian is wanted.
// partial_surface_r_wrt_(rho,sigma)
//	= The partial derivatives of the surface radius with respect to
//	  the (rho,sigma) coordinates.
// g_{xx,xy,xz,yy,yz,zz} = The xyz 3-metric components $g_{ij}$.
// g_{rho_rho,rho_sigma,sigma_sigma} = The (rho,sigma) induced 2-D metric.
//
// Results:
// This function returns the Jacobian of the (rho,sigma) induced 2-D metric.
//
fp patch::rho_sigma_metric(fp r, fp rho, fp sigma,
			   fp partial_surface_r_wrt_rho,
			   fp partial_surface_r_wrt_sigma,
			   fp g_xx, fp g_xy, fp g_xz,
				    fp g_yy, fp g_yz,
					     fp g_zz,
			   fp& g_rho_rho, fp& g_rho_sigma,
					  fp& g_sigma_sigma)
	const
{
fp partial_x_wrt_r, partial_x_wrt_rho, partial_x_wrt_sigma;
fp partial_y_wrt_r, partial_y_wrt_rho, partial_y_wrt_sigma;
fp partial_z_wrt_r, partial_z_wrt_rho, partial_z_wrt_sigma;
partial_xyz_wrt_r_rho_sigma
	(r, rho, sigma,
	 partial_x_wrt_r, partial_x_wrt_rho, partial_x_wrt_sigma,
	 partial_y_wrt_r, partial_y_wrt_rho, partial_y_wrt_sigma,
	 partial_z_wrt_r, partial_z_wrt_rho, partial_z_wrt_sigma);

const fp dx_wrt_rho   = partial_x_wrt_rho
			+ partial_x_wrt_r*partial_surface_r_wrt_rho;
const fp dx_wrt_sigma = partial_x_wrt_sigma
			+ partial_x_wrt_r*partial_surface_r_wrt_sigma;
const fp dy_wrt_rho   = partial_y_wrt_rho
			+ partial_y_wrt_r*partial_surface_r_wrt_rho;
const fp dy_wrt_sigma = partial_y_wrt_sigma
			+ partial_y_wrt_r*partial_surface_r_wrt_sigma;
const fp dz_wrt_rho   = partial_z_wrt_rho
			+ partial_z_wrt_r*partial_surface_r_wrt_rho;
const fp dz_wrt_sigma = partial_z_wrt_sigma
			+ partial_z_wrt_r*partial_surface_r_wrt_sigma;

g_rho_rho     = +     dx_wrt_rho*dx_wrt_rho*g_xx
		+ 2.0*dx_wrt_rho*dy_wrt_rho*g_xy
		+ 2.0*dx_wrt_rho*dz_wrt_rho*g_xz
		+     dy_wrt_rho*dy_wrt_rho*g_yy
		+ 2.0*dy_wrt_rho*dz_wrt_rho*g_yz
		+     dz_wrt_rho*dz_wrt_rho*g_zz;
g_rho_sigma =   +  dx_wrt_rho*dx_wrt_sigma                           *g_xx
		+ (dx_wrt_rho*dy_wrt_sigma + dy_wrt_rho*dx_wrt_sigma)*g_xy
		+ (dx_wrt_rho*dz_wrt_sigma + dz_wrt_rho*dx_wrt_sigma)*g_xz
		+  dy_wrt_rho*dy_wrt_sigma                           *g_yy
		+ (dy_wrt_rho*dz_wrt_sigma + dz_wrt_rho*dy_wrt_sigma)*g_yz
		+  dz_wrt_rho*dz_wrt_sigma                           *g_zz;
g_sigma_sigma = + dx_wrt_sigma*dx_wrt_sigma*g_xx
		+ 2.0*dx_wrt_sigma*dy_wrt_sigma*g_xy
		+ 2.0*dx_wrt_sigma*dz_wrt_sigma*g_xz
		+     dy_wrt_sigma*dy_wrt_sigma*g_yy
		+ 2.0*dy_wrt_sigma*dz_wrt_sigma*g_yz
		+     dz_wrt_sigma*dz_wrt_sigma*g_zz;

return g_rho_rho*g_sigma_sigma - jtutil::pow2(g_rho_sigma);
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function decodes the character-string name of an integration method
// into an  enum integration_method .  See the comments in "patch.hh" on the
// declaration of  enum integration_method  for details on the methods and
// their character-string names.
//
//static
  enum patch::integration_method
    patch::decode_integration_method(const char method_string[])
{
if	(    STRING_EQUAL(method_string, "trapezoid")
	  || STRING_EQUAL(method_string, "trapezoid rule")    )
   then return integration_method__trapezoid;
else if (    STRING_EQUAL(method_string, "Simpson")
	  || STRING_EQUAL(method_string, "Simpson's rule")    )
   then return integration_method__Simpson;
else if (    STRING_EQUAL(method_string, "Simpson (variant)")
	  || STRING_EQUAL(method_string, "Simpson's rule (variant)")    )
   then return integration_method__Simpson_variant;
else if (    STRING_EQUAL(method_string, "automatic choice")    )
   then return integration_method__automatic_choice;
else	error_exit(ERROR_EXIT,
"***** patch::decode_integration_method():\n"
"        unknown method_string=\"%s\"!\n"
,
		   method_string);				/*NOTREACHED*/
}

//******************************************************************************

//
// This function computes an approximation to the arc length of a surface
// over the patch's nominal bounds along the rho direction (i.e. in a
// dsigma=constant plane where dsigma is a multiple of 90 degrees)
//
// Arguments:
// ghosted_radius_gfn = (in) The surface radius.
// g_{xx,xy,xz,yy,yz,zz}_gfn = (in) The xyz 3-metric components.
// method = (in) Selects the integration scheme.
//
fp patch::rho_arc_length(int ghosted_radius_gfn,
			 int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
				       int g_yy_gfn, int g_yz_gfn,
						     int g_zz_gfn,
			 enum integration_method method)
	const
{
fp dsigma;
if	(is_valid_dsigma(  0.0)) then dsigma =   0.0;
else if (is_valid_dsigma( 90.0)) then dsigma =  90.0;
else if (is_valid_dsigma(180.0)) then dsigma = 180.0;
else if (is_valid_dsigma(-90.0)) then dsigma = -90.0;
else	error_exit(PANIC_EXIT,
"***** patch::rho_arc_length(): can't find valid dsigma\n"
"                               which is a multiple of 90 degrees!\n"
"                               %s patch: [min,max]_dsigma()=[%g,%g]\n"
		   ,
		   name(), double(min_dsigma()), double(max_dsigma()));
const fp sigma = sigma_of_dsigma(dsigma);
const int isigma = isigma_of_sigma(sigma);

fp sum = 0.0;

	for (int irho = min_irho() ; irho <= max_irho() ; ++irho)
	{
	const fp rho = rho_of_irho(irho);
	const fp r = ghosted_gridfn(ghosted_radius_gfn, irho,isigma);
	const fp partial_surface_r_wrt_rho
		= partial_rho  (ghosted_radius_gfn, irho,isigma);
	const fp partial_surface_r_wrt_sigma
		= partial_sigma(ghosted_radius_gfn, irho,isigma);

	const fp g_xx = gridfn(g_xx_gfn, irho,isigma);
	const fp g_xy = gridfn(g_xy_gfn, irho,isigma);
	const fp g_xz = gridfn(g_xz_gfn, irho,isigma);
	const fp g_yy = gridfn(g_yy_gfn, irho,isigma);
	const fp g_yz = gridfn(g_yz_gfn, irho,isigma);
	const fp g_zz = gridfn(g_zz_gfn, irho,isigma);

	fp g_rho_rho, g_rho_sigma, g_sigma_sigma;
	rho_sigma_metric(r, rho, sigma,
			 partial_surface_r_wrt_rho,
			 partial_surface_r_wrt_sigma,
			 g_xx, g_xy, g_xz,
			       g_yy, g_yz,
				     g_zz,
			 g_rho_rho, g_rho_sigma,
				    g_sigma_sigma);

	const fp coeff = integration_coeff(method,
					   max_irho()-min_irho(),
					   irho      -min_irho());

	sum += coeff * sqrt(g_rho_rho);
	}

return delta_rho() * sum;
}

//******************************************************************************

//
// This function computes an approximation to the arc length of a surface
// over the patch's nominal bounds along the sigma direction (i.e. in a
// drho=constant plane where drho is a multiple of 90 degrees)
//
// Arguments:
// ghosted_radius_gfn = (in) The surface radius.
// g_{xx,xy,xz,yy,yz,zz}_gfn = (in) The xyz 3-metric components.
// method = (in) Selects the integration scheme.
//
fp patch::sigma_arc_length(int ghosted_radius_gfn,
			   int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
					 int g_yy_gfn, int g_yz_gfn,
						       int g_zz_gfn,
			   enum integration_method method)
	const
{
fp drho;
if	(is_valid_drho(  0.0)) then drho =   0.0;
else if (is_valid_drho( 90.0)) then drho =  90.0;
else if (is_valid_drho(180.0)) then drho = 180.0;
else if (is_valid_drho(-90.0)) then drho = -90.0;
else	error_exit(PANIC_EXIT,
"***** patch::sigma_arc_length(): can't find valid drho\n"
"                                 which is a multiple of 90 degrees!\n"
"                                 %s patch: [min,max]_drho()=[%g,%g]\n"
		   ,
		   name(), double(min_drho()), double(max_drho()));
const fp rho = rho_of_drho(drho);
const int irho = irho_of_rho(rho);

fp sum = 0.0;

	for (int isigma = min_isigma() ; isigma <= max_isigma() ; ++isigma)
	{
	const fp sigma = sigma_of_isigma(isigma);
	const fp r = ghosted_gridfn(ghosted_radius_gfn, irho,isigma);
	const fp partial_surface_r_wrt_rho
		= partial_rho  (ghosted_radius_gfn, irho,isigma);
	const fp partial_surface_r_wrt_sigma
		= partial_sigma(ghosted_radius_gfn, irho,isigma);

	const fp g_xx = gridfn(g_xx_gfn, irho,isigma);
	const fp g_xy = gridfn(g_xy_gfn, irho,isigma);
	const fp g_xz = gridfn(g_xz_gfn, irho,isigma);
	const fp g_yy = gridfn(g_yy_gfn, irho,isigma);
	const fp g_yz = gridfn(g_yz_gfn, irho,isigma);
	const fp g_zz = gridfn(g_zz_gfn, irho,isigma);

	fp g_rho_rho, g_rho_sigma, g_sigma_sigma;
	rho_sigma_metric(r, rho, sigma,
			 partial_surface_r_wrt_rho,
			 partial_surface_r_wrt_sigma,
			 g_xx, g_xy, g_xz,
			       g_yy, g_yz,
				     g_zz,
			 g_rho_rho, g_rho_sigma,
				    g_sigma_sigma);

	const fp coeff = integration_coeff(method,
					   max_isigma()-min_isigma(),
					   isigma      -min_isigma());

	sum += coeff * sqrt(g_sigma_sigma);
	}

return delta_sigma() * sum;
}

//******************************************************************************

//
// This function computes the arc length of a surface in the specified
// plane ("xz" or "yz") over the patch's nominal bounds.
//
// Arguments:
// plane[] = (in) "xz" or "yz" to specify the plane.
// ghosted_radius_gfn = (in) The surface radius.
// g_{xx,xy,xz,yy,yz,zz}_gfn = (in) The xyz 3-metric components.
// method = (in) Selects the integration scheme.
//
fp z_patch::plane_arc_length(const char plane[],
			     int ghosted_radius_gfn,
			     int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
					   int g_yy_gfn, int g_yz_gfn,
							 int g_zz_gfn,
			     enum integration_method method)
	const
{
if	((plane[0] == 'x') && (plane[1] == 'z'))
   then // xz-plane = rotation about y = nu arc = sigma sigma
	return sigma_arc_length(ghosted_radius_gfn,
				g_xx_gfn, g_xy_gfn, g_xz_gfn,
					  g_yy_gfn, g_yz_gfn,
						    g_zz_gfn,
				method);
else if ((plane[0] == 'y') && (plane[1] == 'z'))
   then // yz-plane = rotation about x = mu arc = rho arc
	return rho_arc_length(ghosted_radius_gfn,
			      g_xx_gfn, g_xy_gfn, g_xz_gfn,
					g_yy_gfn, g_yz_gfn,
						  g_zz_gfn,
			      method);
else	error_exit(ERROR_EXIT,
"***** z_patch::plane_arc_length(): %s patch, plane=\"%s\", but\n"
"                                   this patch doesn't contain that plane!\n"
,
		   name(), plane);				/*NOTREACHED*/
}

//******************************************************************************

//
// This function computes the arc length of a surface in the specified
// plane ("xy" or "xz") over the patch's nominal bounds.
//
// Arguments:
// plane[] = (in) "xy" or "xz" to specify the plane.
// ghosted_radius_gfn = (in) The surface radius.
// g_{xx,xy,xz,yy,yz,zz}_gfn = (in) The xyz 3-metric components.
// method = (in) Selects the integration scheme.
//
fp x_patch::plane_arc_length(const char plane[],
			     int ghosted_radius_gfn,
			     int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
					   int g_yy_gfn, int g_yz_gfn,
							 int g_zz_gfn,
			     enum integration_method method)
	const
{
if	((plane[0] == 'x') && (plane[1] == 'y'))
   then // xy-plane = rotation about z = phi arc = sigma arc
	return sigma_arc_length(ghosted_radius_gfn,
				g_xx_gfn, g_xy_gfn, g_xz_gfn,
					  g_yy_gfn, g_yz_gfn,
						    g_zz_gfn,
				method);
else if ((plane[0] == 'x') && (plane[1] == 'z'))
   then // xz-plane = rotation about y = nu arc = rho arc
	return rho_arc_length(ghosted_radius_gfn,
			      g_xx_gfn, g_xy_gfn, g_xz_gfn,
					g_yy_gfn, g_yz_gfn,
						  g_zz_gfn,
			      method);
else	error_exit(ERROR_EXIT,
"***** x_patch::plane_arc_length(): %s patch, plane=\"%s\", but\n"
"                                   this patch doesn't contain that plane!\n"
,
		   name(), plane);				/*NOTREACHED*/
}

//******************************************************************************

//
// This function computes the arc length of a surface in the specified
// plane ("xy" or "yz") over the patch's nominal bounds.
//
// Arguments:
// plane[] = (in) "xy" or "yz" to specify the plane.
// ghosted_radius_gfn = (in) The surface radius.
// g_{xx,xy,xz,yy,yz,zz}_gfn = (in) The xyz 3-metric components.
// method = (in) Selects the integration scheme.
//
fp y_patch::plane_arc_length(const char plane[],
			     int ghosted_radius_gfn,
			     int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
					   int g_yy_gfn, int g_yz_gfn,
							 int g_zz_gfn,
			     enum integration_method method)
	const
{
if	((plane[0] == 'x') && (plane[1] == 'y'))
   then // xy-plane = rotation about z = phi arc = sigma arc
	return sigma_arc_length(ghosted_radius_gfn,
				g_xx_gfn, g_xy_gfn, g_xz_gfn,
					  g_yy_gfn, g_yz_gfn,
						    g_zz_gfn,
				method);
else if ((plane[0] == 'y') && (plane[1] == 'z'))
   then // yz-plane = rotation about x = mu arc = rho arc
	return rho_arc_length(ghosted_radius_gfn,
			      g_xx_gfn, g_xy_gfn, g_xz_gfn,
					g_yy_gfn, g_yz_gfn,
						  g_zz_gfn,
			      method);
else	error_exit(ERROR_EXIT,
"***** y_patch::plane_arc_length(): %s patch, plane=\"%s\", but\n"
"                                   this patch doesn't contain that plane!\n"
,
		   name(), plane);				/*NOTREACHED*/
}

//******************************************************************************

//
// This function computes an approximation to the (surface) integral of
// a gridfn over the patch's nominal area,
//	$\int f(\rho,\sigma) \, dA$
//		= \int f(\rho,\sigma) \sqrt{|J|} \, d\rho \, d\sigma$
// where $J$ is the Jacobian of $(x,y,z)$ with respect to $(rho,sigma).
//
// Arguments:
// unknown_src_gfn = (in) The gridfn to be integrated.  This may be
//			  either nominal-grid or ghosted-grid; n.b. in
//			  the latter case the integral is still done only
//			  over the patch's nominal area.
// ghosted_radius_gfn = (in) The surface radius.
// g_{xx,xy,xz,yy,yz,zz}_gfn = (in) The xyz 3-metric components.
// method = (in) Selects the integration scheme.
//
fp patch::integrate_gridfn(int unknown_src_gfn,
			   int ghosted_radius_gfn,
			   int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
					 int g_yy_gfn, int g_yz_gfn,
						       int g_zz_gfn,
			   enum integration_method method)
	const
{
const bool src_is_ghosted = is_valid_ghosted_gfn(unknown_src_gfn);

fp sum = 0.0;
	for (int irho = min_irho() ; irho <= max_irho() ; ++irho)
	{
	for (int isigma = min_isigma() ; isigma <= max_isigma() ; ++isigma)
	{
	const fp fn = unknown_gridfn(src_is_ghosted,
				     unknown_src_gfn, irho,isigma);

	const fp rho   = rho_of_irho    (irho);
	const fp sigma = sigma_of_isigma(isigma);
	const fp r = ghosted_gridfn(ghosted_radius_gfn, irho,isigma);
	const fp partial_surface_r_wrt_rho
		= partial_rho  (ghosted_radius_gfn, irho,isigma);
	const fp partial_surface_r_wrt_sigma
		= partial_sigma(ghosted_radius_gfn, irho,isigma);

	const fp g_xx = gridfn(g_xx_gfn, irho,isigma);
	const fp g_xy = gridfn(g_xy_gfn, irho,isigma);
	const fp g_xz = gridfn(g_xz_gfn, irho,isigma);
	const fp g_yy = gridfn(g_yy_gfn, irho,isigma);
	const fp g_yz = gridfn(g_yz_gfn, irho,isigma);
	const fp g_zz = gridfn(g_zz_gfn, irho,isigma);

	fp g_rho_rho, g_rho_sigma, g_sigma_sigma;
	const fp Jac = rho_sigma_metric(r, rho, sigma,
					partial_surface_r_wrt_rho,
					partial_surface_r_wrt_sigma,
					g_xx, g_xy, g_xz,
					      g_yy, g_yz,
						    g_zz,
					g_rho_rho, g_rho_sigma,
						   g_sigma_sigma);

	const fp coeff_rho   = integration_coeff(method,
						 max_irho()-min_irho(),
						 irho      -min_irho());
	const fp coeff_sigma = integration_coeff(method,
						 max_isigma()-min_isigma(),
						 isigma      -min_isigma());

	sum += coeff_rho*coeff_sigma * fn * sqrt(jtutil::abs(Jac));
	}
	}

return delta_rho() * delta_sigma() * sum;
}

fp patch::integrate_gridpoint(int unknown_src_gfn,
                              int ghosted_radius_gfn,
                              int g_xx_gfn, int g_xy_gfn, int g_xz_gfn,
					    int g_yy_gfn, int g_yz_gfn,
                                                          int g_zz_gfn,
                              enum integration_method method,
                              int irho, int isigma)
	const
{
const bool src_is_ghosted = is_valid_ghosted_gfn(unknown_src_gfn);

const fp fn = unknown_gridfn(src_is_ghosted,
			     unknown_src_gfn, irho,isigma);

const fp rho   = rho_of_irho    (irho);
const fp sigma = sigma_of_isigma(isigma);
const fp r = ghosted_gridfn(ghosted_radius_gfn, irho,isigma);
const fp partial_surface_r_wrt_rho
	= partial_rho  (ghosted_radius_gfn, irho,isigma);
const fp partial_surface_r_wrt_sigma
	= partial_sigma(ghosted_radius_gfn, irho,isigma);

const fp g_xx = gridfn(g_xx_gfn, irho,isigma);
const fp g_xy = gridfn(g_xy_gfn, irho,isigma);
const fp g_xz = gridfn(g_xz_gfn, irho,isigma);
const fp g_yy = gridfn(g_yy_gfn, irho,isigma);
const fp g_yz = gridfn(g_yz_gfn, irho,isigma);
const fp g_zz = gridfn(g_zz_gfn, irho,isigma);

fp g_rho_rho, g_rho_sigma, g_sigma_sigma;
const fp Jac = rho_sigma_metric(r, rho, sigma,
				partial_surface_r_wrt_rho,
				partial_surface_r_wrt_sigma,
				g_xx, g_xy, g_xz,
				      g_yy, g_yz,
					    g_zz,
				g_rho_rho, g_rho_sigma,
					   g_sigma_sigma);

const fp coeff_rho   = integration_coeff(method,
					 max_irho()-min_irho(),
					 irho      -min_irho());
const fp coeff_sigma = integration_coeff(method,
					 max_isigma()-min_isigma(),
					 isigma      -min_isigma());

const fp val = coeff_rho*coeff_sigma * fn * sqrt(jtutil::abs(Jac));

return delta_rho() * delta_sigma() * val;
}

//******************************************************************************

//
// This function computes the integration coefficients for
//  integrate_gridfn() .  That is, if we write
//	$\int_{x_0}^{x_N} f(x) \, dx
//		\approx \Delta x \, \sum_{i=0}^N c_i f(x_i)$
// then this function computes $c_i$.
//
// For method == integration_method__automatic_choice the choices are
//	N=1		trapezoid
//	N=2		Simpson
//	N=3		trapezoid
//	N=4		Simpson
//	N=5		trapezoid
//	N=6		Simpson
//	N=7 and up	Simpson variant
//
// Arguments:
// method = Specifies the integration method.
// N = The number of integration *intervals*.  (The number of integration
//     *points* is N+1.)
// i = Specifies the point at which the coefficient is desired.
//
//static
  fp patch::integration_coeff(enum integration_method method, int N, int i)
{
TBOX_ASSERT(i >= 0);
TBOX_ASSERT(i <= N);

if (method == integration_method__automatic_choice)
   then {
	if	(N >= 7)
	   then method = integration_method__Simpson_variant;
	else if ((N % 2) == 0)
	   then method = integration_method__Simpson;
	else	method = integration_method__trapezoid;
	}

switch	(method)
	{
case integration_method__trapezoid:
	if ((i == 0) || (i == N))
	   then return 0.5;
	   else return 1.0;

case integration_method__Simpson:
	if ((N % 2) != 0)
	   then error_exit(ERROR_EXIT,
"***** patch::integration_coeff():\n"
"        Simpson's rule requires N to be even, but N=%d!\n",
			   N);					/*NOTREACHED*/
	if	((i == 0) || (i == N))
	   then return 1.0/3.0;
	else if ((i % 2) == 0)
	   then return 2.0/3.0;
	else	return 4.0/3.0;

case integration_method__Simpson_variant:
	if (N < 7)
	   then error_exit(ERROR_EXIT,
"***** patch::integration_coeff():\n"
"        Simpson's rule (variant) requires N >= 7, but N=%d!\n",
			   N);					/*NOTREACHED*/
	if	((i == 0) || (i == N))
	   then return 17.0/48.0;
	else if ((i == 1) || (i == N-1))
	   then return 59.0/48.0;
	else if ((i == 2) || (i == N-2))
	   then return 43.0/48.0;
	else if ((i == 3) || (i == N-3))
	   then return 49.0/48.0;
	else	return 1.0;

default:
	error_exit(ERROR_EXIT,
"***** patch::integration_coeff(): unknown method=(int)%d!\n"
"                                  (this should never happen!)\n"
,
		   int(method));				/*NOTREACHED*/
	}
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function returns a reference to the ghost zone on a specified
// edge, after first TBOX_ASSERT()ing that the edge belongs to this patch.
//
// N.b. This function can't be inline in "patch.hh" because it needs
//	member functions of class patch_edge, which comes after class patch
//	in our #include order.
//
ghost_zone& patch::ghost_zone_on_edge(const patch_edge& e)
	const
{
TBOX_ASSERT(e.my_patch() == *this);
return minmax_ang_ghost_zone(e.is_min(), e.is_rho());
}

//******************************************************************************

//
// This function determines which of the two adjacent ghost zones meeting
// at a specified corner, contains a specified point.  If the point isn't
// in either ghost zone, an error_exit() is done.  If the point is in both
// ghost zones, it's arbitrary which one will be chosen.
//
// Arguments:
// {rho,sigma}_is_min = Specify the corner (and implicitly the ghost zones).
// irho,isigma = Specify the point.
//
// Results:
// This function returns (a reference to) the desired ghost zone.
ghost_zone& patch::corner_ghost_zone_containing_point
	(bool rho_is_min, bool sigma_is_min,
	 int irho, int isigma)
	const
{
ghost_zone&   rho_gz =   minmax_rho_ghost_zone(  rho_is_min);
ghost_zone& sigma_gz = minmax_sigma_ghost_zone(sigma_is_min);

const patch_edge&   rho_edge =   rho_gz.my_edge();
const patch_edge& sigma_edge = sigma_gz.my_edge();

const int   rho_iperp =   rho_edge.iperp_of_irho_isigma(irho, isigma);
const int   rho_ipar  =   rho_edge. ipar_of_irho_isigma(irho, isigma);
const int sigma_iperp = sigma_edge.iperp_of_irho_isigma(irho, isigma);
const int sigma_ipar  = sigma_edge. ipar_of_irho_isigma(irho, isigma);

const bool is_in_rho_ghost_zone
	=   rho_gz.is_in_ghost_zone(  rho_iperp,   rho_ipar);
const bool is_in_sigma_ghost_zone
	= sigma_gz.is_in_ghost_zone(sigma_iperp, sigma_ipar);

// check that point is in at least one ghost zone
if (!is_in_rho_ghost_zone && !is_in_sigma_ghost_zone)
   then error_exit(ERROR_EXIT,
"***** patch::corner_ghost_zone_containing_point():\n"
"        neither ghost zone contains point (this should never happen)!\n"
"        patch=%s rho_is_min=(int)%d sigma_is_min=(int)%d\n"
"        irho=%d isigma=%d\n"
,
	   name(), int(rho_is_min), int(sigma_is_min),
	   irho, isigma);					/*NOTREACHED*/

return is_in_rho_ghost_zone ? rho_gz : sigma_gz;
}

//******************************************************************************

//
// This function determines which ghost zone contains a specified
// noncorner point.
//
// If the point isn't in any ghost zone of this patch, or if the point
// is in the corner of a ghost zone, an error_exit() is done.
//
// Arguments:
// irho,isigma = Specify the point.
//
// Results:
// This function returns (a reference to) the desired ghost zone.
ghost_zone& patch::ghost_zone_containing_noncorner_point(int irho, int isigma)
	const
{
	// n.b. these loops must use _int_ variables for the loop
	//      to terminate!
	for (int want_min = false ; want_min <= true ; ++want_min)
	{
	for (int want_rho = false ; want_rho <= true ; ++want_rho)
	{
	const patch_edge& e = minmax_ang_patch_edge(want_min, want_rho);
	const int iperp = e.iperp_of_irho_isigma(irho, isigma);
	const int ipar  = e.ipar_of_irho_isigma (irho, isigma);

	ghost_zone& gz = minmax_ang_ghost_zone(want_min, want_rho);
	if ( gz.is_in_ghost_zone(iperp, ipar)
	     && gz.my_edge().ipar_is_in_noncorner(ipar) )
	   then return gz;
	}
	}

error_exit(ERROR_EXIT,
"***** patch::ghost_zone_containing_noncorner_point():\n"
"        no ghost zone contains point (this should never happen)!\n"
"        patch=%s irho=%d isigma=%d\n"
,
	   name(), irho, isigma);				/*NOTREACHED*/
}

//******************************************************************************

//
// This function TBOX_ASSERT()s that a specified ghost zone of this patch
// hasn't already been set up, then constructs it as a mirror-symmetry
// ghost zone and properly links this to/from the patch.
//
void patch::create_mirror_symmetry_ghost_zone(const patch_edge& my_edge)
{
// make sure we belong to the right patch
TBOX_ASSERT(my_edge.my_patch() == *this);

symmetry_ghost_zone *temp = new symmetry_ghost_zone(my_edge);
set_ghost_zone(my_edge, temp);
}

//******************************************************************************

//
// This function TBOX_ASSERT()s that a specified ghost zone of this patch
// hasn't already been set up, then creates it as a periodic-symmetry
// ghost zone and properly links this to/from the patch.
//
void patch::create_periodic_symmetry_ghost_zone
	(const patch_edge& my_edge, const patch_edge& other_edge,
	 bool is_ipar_map_plus)
{
// make sure we belong to the right patch
TBOX_ASSERT(my_edge.my_patch() == *this);

int my_sample_ipar = my_edge.min_ipar_without_corners();
int other_sample_ipar = is_ipar_map_plus
			? other_edge.min_ipar_without_corners()
			: other_edge.max_ipar_without_corners();

symmetry_ghost_zone *temp
	= new symmetry_ghost_zone(my_edge,        other_edge,
				  my_sample_ipar, other_sample_ipar,
				  is_ipar_map_plus);
set_ghost_zone(my_edge, temp);
}

//******************************************************************************

//
// This function TBOX_ASSERT()s that a specified ghost zone of this patch
// hasn't already been set up, then creates it as an interpatch ghost
// zone (with lots of NULL pointers for info we can't compute yet)
// and properly links this to/from the patch.
//
void patch::create_interpatch_ghost_zone
	(const patch_edge& my_edge, const patch_edge& other_edge,
	 int patch_overlap_width)
{
// make sure we belong to the right patch
TBOX_ASSERT(my_edge.my_patch() == *this);

interpatch_ghost_zone *temp
	= new interpatch_ghost_zone(my_edge, other_edge,
				    patch_overlap_width);
set_ghost_zone(my_edge, temp);
}

//******************************************************************************

//
// This is a helper function for  setup_*_ghost_zone().  This function
// TBOX_ASSERT()s that one of the ghost zone pointers (which one is selected
// by  edge ) is NULL, then stores a value in it.
//
void patch::set_ghost_zone(const patch_edge& edge, ghost_zone* gzp)
{
ghost_zone*& ghost_zone_ptr_to_set
	= edge.is_min()
	  ? (edge.is_rho() ? min_rho_ghost_zone_ : min_sigma_ghost_zone_)
	  : (edge.is_rho() ? max_rho_ghost_zone_ : max_sigma_ghost_zone_);

TBOX_ASSERT(ghost_zone_ptr_to_set == NULL);
ghost_zone_ptr_to_set = gzp;
}

//******************************************************************************

//
// This function finds which patch edge is adjacent to a neighboring
// patch q, or does an error_exit() if q isn't actually a neighboring patch.
// The computation is done using only (rho,sigma) coordinate sets and
// min/max dang bounds ==> it's ok to use this function in setting up
// interpatch ghost zones.
//
// Arguments:
// q = The (supposedly) neighboring patch.
// patch_overlap_width = The number of grid points these patches overlap.
//		      If this is nonzero, then these patches must have the
//		      same grid spacing in the perpendicular direction.
//
const patch_edge& patch::edge_adjacent_to_patch(const patch& q,
						int patch_overlap_width /* = 0 */)
	const
{
const patch& p = *this;

// which (rho,sigma) coordinate do the patches have in common?
// ... this is the perp coordinate for the border
const local_coords::coords_set common_coord_set
	= p.coords_set_rho_sigma() & q.coords_set_rho_sigma();

// is this coordinate rho or sigma in each patch?
const bool common_is_p_rho   = (common_coord_set == p.coords_set_rho  ());
const bool common_is_p_sigma = (common_coord_set == p.coords_set_sigma());
if ((common_is_p_rho ^ common_is_p_sigma) != 0x1)
   then error_exit(ERROR_EXIT,
"***** patch::edge_adjacent_to_patch():\n"
"        common coordinate isn't exactly one of p.{rho,sigma}!\n"
"        p.name()=\"%s\" q.name()=\"%s\"\n"
"        common_coord_set=%s\n"
"        common_is_p_rho=%d common_is_p_sigma=%d\n"
,
		   p.name(), q.name(),
		   local_coords::name_of_coords_set(common_coord_set),
		   int(common_is_p_rho), int(common_is_p_sigma));
								/*NOTREACHED*/
const bool common_is_q_rho   = (common_coord_set == q.coords_set_rho  ());
const bool common_is_q_sigma = (common_coord_set == q.coords_set_sigma());
if ((common_is_q_rho ^ common_is_q_sigma) != 0x1)
   then error_exit(ERROR_EXIT,
"***** patch::edge_adjacent_to_patch():\n"
"        common coordinate isn't exactly one of q.{rho,sigma}!\n"
"        p.name()=\"%s\" q.name()=\"%s\"\n"
"        common_coord_set=%s\n"
"        common_is_q_rho=%d common_is_q_sigma=%d\n"
,
		   p.name(), q.name(),
		   local_coords::name_of_coords_set(common_coord_set),
		   int(common_is_q_rho), int(common_is_q_sigma));
								/*NOTREACHED*/

// how much do the patches overlap?
// ... eg patch_overlap_width = 3 would be
//	p   p   p   p   p
//		q   q   q   q   q
//     so the overlap would be (patch_overlap_width-1) * delta = 2 * delta
if ( (patch_overlap_width-1 != 0)
     && jtutil::fuzzy<fp>::NE(p.delta_dang(common_is_p_rho),
			      q.delta_dang(common_is_q_rho)) )
   then error_exit(ERROR_EXIT,
"***** patch::edge_adjacent_to_patch():\n"
"        patch_overlap_width != 0 must have same perp grid spacing in both patches!\n"
"        p.name()=\"%s\" q.name()=\"%s\"\n"
"        common_coord_set=%s\n"
"        common_is_p_rho=%d common_is_q_rho=%d\n"
"        p.delta_dang(common_is_p_rho)=%g\n"
"        q.delta_dang(common_is_q_rho)=%g\n"
,
		   p.name(), q.name(),
		   local_coords::name_of_coords_set(common_coord_set),
		   int(common_is_p_rho), int(common_is_q_rho),
		   double(p.delta_dang(common_is_p_rho)),
		   double(q.delta_dang(common_is_q_rho)));	/*NOTREACHED*/


const fp doverlap = fp(patch_overlap_width-1) * p.delta_dang(common_is_p_rho);

// where is the common boundary relative to the min/max sides of each patch?
const bool common_is_p_min_q_max
    = local_coords::fuzzy_EQ_dang(p.min_dang(common_is_p_rho),
				  q.max_dang(common_is_q_rho) - doverlap);
const bool common_is_p_max_q_min
    = local_coords::fuzzy_EQ_dang(p.max_dang(common_is_p_rho),
				  q.min_dang(common_is_q_rho) + doverlap);
if ((common_is_p_min_q_max ^ common_is_p_max_q_min) != 0x1)
   then error_exit(ERROR_EXIT,
"***** patch::edge_adjacent_to_patch():\n"
"        common coordinate isn't exactly one of {pmax/qmin, pmin/qmax}!\n"
"        p.name()=\"%s\" q.name()=\"%s\"\n"
"        common_coord_set=%s\n"
"        common_is_p_rho=%d common_is_q_rho=%d\n"
"        p.delta_dang(common_is_p_rho)=%g\n"
"        q.delta_dang(common_is_q_rho)=%g\n"
"        patch_overlap_width=%d doverlap=%g\n"
"        common_is_p_min_q_max=%d common_is_p_max_q_min=%d\n"
,
		   p.name(), q.name(),
		   local_coords::name_of_coords_set(common_coord_set),
		   int(common_is_p_rho), int(common_is_q_rho),
		   double(p.delta_dang(common_is_p_rho)),
		   double(q.delta_dang(common_is_q_rho)),
		   patch_overlap_width, double(doverlap),
		   int(common_is_p_min_q_max), int(common_is_p_max_q_min));
								/*NOTREACHED*/

return p.minmax_ang_patch_edge(common_is_p_min_q_max, common_is_p_rho);
}

//******************************************************************************

//
// This function verifies (via TBOX_ASSERT()) that all ghost zones of this
// patch have been fully set up.
//
void patch::assert_all_ghost_zones_fully_setup() const
{
TBOX_ASSERT(min_rho_ghost_zone_ != NULL);
TBOX_ASSERT(max_rho_ghost_zone_ != NULL);
TBOX_ASSERT(min_sigma_ghost_zone_ != NULL);
TBOX_ASSERT(max_sigma_ghost_zone_ != NULL);

// these calls are no-ops for non-interpatch ghost zones
min_rho_ghost_zone().assert_fully_setup();
max_rho_ghost_zone().assert_fully_setup();
min_sigma_ghost_zone().assert_fully_setup();
max_sigma_ghost_zone().assert_fully_setup();
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

	  }	// namespace AHFinderDirect
