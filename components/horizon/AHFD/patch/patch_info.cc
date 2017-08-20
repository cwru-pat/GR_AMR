// patch_info.cc -- POD struct of minimal info varying from one patch to another
// $Header$

//
// patch_info::grid_array_pars
// patch_info::grid_pars
/// patch_info::verify_grid_spacing_ok
//

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "../cctk.h"

#include "../jtutil/util.hh"
#include "../jtutil/array.hh"
#include "../jtutil/cpm_map.hh"
#include "../jtutil/linear_map.hh"

#include "coords.hh"
#include "grid.hh"
#include "patch_info.hh"

// all the code in this file is inside this namespace
namespace AHFD
	  {

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function computes, and returns a reference to, a
//  struct grid_arrays::grid_array_pars  from the info in a
//  struct patch_info  and the additional information in the arguments.
//
// The result refers to an internal static buffer in this function; the
// usual caveats about lifetimes/overwriting apply.
//
// Arguments:
// ghost_zone_width = Width in grid points of all ghost zones.
// patch_extend_width = Number of grid points to extend each patch past
//		     "just touching" so as to overlap neighboring patches.
//		     Thus patches overlap by
//			patch_overlap_width = 2*patch_extend_width + 1
//		     grid points.  For example, with patch_extend_width == 2,
//		     here are the grid points of two neighboring patches:
//			x   x   x   x   x   X   X
//                                      |
//			        O   O   o   o   o   o   o
//		     Here | marks the "just touching" boundary,
//		     x and o the grid points before this extension,
//		     and X and O the extra grid points added by this
//		     extension.
// N_zones_per_right_angle = This sets the grid spacing (same in both
//			     directions) to 90.0 / N_zones_per_right_angle.
//			     It's a fatal error (error_exit()) if this
//			     doesn't evenly divide the grid sizes in both
//			     directions.
//
const grid_arrays::grid_array_pars&
  patch_info::grid_array_pars(int ghost_zone_width, int patch_extend_width,
			      int N_zones_per_right_angle)
	const
{
static
  struct grid_arrays::grid_array_pars grid_array_pars_buffer;

//
// the values of min_(irho,isigma) are actually arbitrary, but for
// debugging convenience it's handy to have (irho,isigma) ranges map
// one-to-one with (rho,sigma) ranges across all patches; the assignments
// here have this property
//
const fp delta_drho_dsigma = 90.0 / fp(N_zones_per_right_angle);
grid_array_pars_buffer.min_irho
	= jtutil::round<fp>::to_integer(min_drho  /delta_drho_dsigma);
grid_array_pars_buffer.min_isigma
	= jtutil::round<fp>::to_integer(min_dsigma/delta_drho_dsigma);

verify_grid_spacing_ok(N_zones_per_right_angle);
const int N_irho_zones
	= jtutil::round<fp>::to_integer(
		   fp(N_zones_per_right_angle) * (max_drho  -min_drho  ) / 90.0
				       );
const int N_isigma_zones
	= jtutil::round<fp>::to_integer(
		   fp(N_zones_per_right_angle) * (max_dsigma-min_dsigma) / 90.0
				       );

grid_array_pars_buffer.max_irho
	= grid_array_pars_buffer.min_irho   + N_irho_zones;
grid_array_pars_buffer.max_isigma
	= grid_array_pars_buffer.min_isigma + N_isigma_zones;

grid_array_pars_buffer.min_irho   -= patch_extend_width;
grid_array_pars_buffer.min_isigma -= patch_extend_width;
grid_array_pars_buffer.max_irho   += patch_extend_width;
grid_array_pars_buffer.max_isigma += patch_extend_width;

grid_array_pars_buffer.min_rho_ghost_zone_width = ghost_zone_width;
grid_array_pars_buffer.max_rho_ghost_zone_width = ghost_zone_width;
grid_array_pars_buffer.min_sigma_ghost_zone_width = ghost_zone_width;
grid_array_pars_buffer.max_sigma_ghost_zone_width = ghost_zone_width;

return grid_array_pars_buffer;
}

//******************************************************************************
//
//
// This function computes, and returns a reference to, a
//  struct grid_arrays::grid_pars  from the info in a  struct patch_info
// and the additional information in the arguments.
//
// The result refers to an internal static buffer in this function; the
// usual caveats about lifetimes/overwriting apply.
//
// Arguments:
// patch_extend_width = Number of grid points to extend each patch past
//		     "just touching" so as to overlap neighboring patches.
//		     Thus patches overlap by  2*patch_extend_width + 1  grid
//		     points.  For example, with patch_extend_width == 2, here
//		     are the grid points of two neighboring patches:
//			x   x   x   x   x   X   X
//                                      |
//			        O   O   o   o   o   o   o
//		     Here | marks the "just touching" boundary,
//		     x and o the grid points before this extension,
//		     and X and O the extra grid points added by this
//		     extension.
// N_zones_per_right_angle = This sets the grid spacing (same in both
//			     directions) to 90.0 / N_zones_per_right_angle.
//			     It's a fatal error (error_exit()) if this
//			     doesn't evenly divide the grid sizes in both
//			     directions.
//
const grid::grid_pars& patch_info::grid_pars(int patch_extend_width,
					     int N_zones_per_right_angle)
	const
{
static
  struct grid::grid_pars grid_pars_buffer;

verify_grid_spacing_ok(N_zones_per_right_angle);
const fp delta_drho_dsigma = 90.0 / fp(N_zones_per_right_angle);
const fp extend_drho_dsigma = fp(patch_extend_width) * delta_drho_dsigma;

grid_pars_buffer.  min_drho   = min_drho   - extend_drho_dsigma;
grid_pars_buffer.delta_drho   = delta_drho_dsigma;
grid_pars_buffer.  max_drho   = max_drho   + extend_drho_dsigma;
grid_pars_buffer.  min_dsigma = min_dsigma - extend_drho_dsigma;
grid_pars_buffer.delta_dsigma = delta_drho_dsigma;
grid_pars_buffer.  max_dsigma = max_dsigma + extend_drho_dsigma;

return grid_pars_buffer;
}

//******************************************************************************

//
// This function verifies that the grid spacing evenly divides the
// grid sizes in both directions, and does an  error_exit()  if not.
//
// Arguments:
// N_zones_per_right_angle = This sets the grid spacing (same in both
//			     directions) to 90.0 / N_zones_per_right_angle.
//
void patch_info::verify_grid_spacing_ok(int N_zones_per_right_angle)
	const
{
const fp N_irho_zones_fp
	= fp(N_zones_per_right_angle) * (max_drho  -min_drho  ) / 90.0;
const fp N_isigma_zones_fp
	= fp(N_zones_per_right_angle) * (max_dsigma-min_dsigma) / 90.0;

if (! (    jtutil::fuzzy<fp>::is_integer(N_irho_zones_fp)
	&& jtutil::fuzzy<fp>::is_integer(N_isigma_zones_fp)    ) )
   then error_exit(ERROR_EXIT,
"***** patch_info::verify_grid_spacing_ok():\n"
"        N_zones_per_right_angle=%d gives grid spacing which\n"
"        doesn't evenly divide grid sizes!\n"
"        [min,max]_drho=[%g,%g] [min,max]_dsigma=[%g,%g]\n"
"        ==> N_irho_zones_fp=%g N_isigma_zones_fp=%g\n"
		   ,
		   N_zones_per_right_angle,
		   double(min_drho), double(max_drho),
		   double(min_dsigma), double(max_dsigma),
		   double(N_irho_zones_fp), double(N_isigma_zones_fp));
								/*NOTREACHED*/
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

	  }	// namespace AHFinderDirect
