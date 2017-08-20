#ifndef AHFD_PATCH_PATCH_INFO_H
#define AHFD_PATCH_PATCH_INFO_H
// patch_info.hh -- POD struct of minimal info varying from one patch to another
// $Header$

//
// prerequisites:
//    <stdio.h>
//    <assert.h>
//    <math.h>
//    "../jtutil/util.hh"
//    "../jtutil//array.hh"
//    "../jtutil/linear_map.hh"
//    "coords.hh"
//    "grid.hh"
//

// everything in this file is inside this namespace
namespace AHFD
	  {

//*****************************************************************************

//
// This (POD, and hence static-initializable) struct gives a minimal
// set of information which varies from one patch to another.
//
// The member functions allow computing all the grid:: constructor
// arguments; with these in hand it's fairly easy to construct the
// patch itself.  This scheme doesn't allow the most general possible
// type of patch (eg it constrains all ghost zones to have the same width,
// and it requires the grid spacing to evenly divide 90 degrees), but
// it does cover all the cases that seem to come up in practice.
//
// Arguments for member functions:
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
struct	patch_info
	{
	const char* name;
	bool is_plus;
	char ctype;
	fp min_drho, max_drho;
	fp min_dsigma, max_dsigma;

	// compute and return reference to  struct grid_arrays::grid_array_pars
	// ... result refers to internal static buffer;
	//     the usual caveats about lifetimes/overwriting apply
	const grid_arrays::grid_array_pars&
	  grid_array_pars(int ghost_zone_width, int patch_extend_width,
			  int N_zones_per_right_angle)
		const;

	// compute and return reference to  struct grid::grid_pars
	// ... result refers to internal static buffer;
	//     the usual caveats about lifetimes/overwriting apply
	const grid::grid_pars& grid_pars(int patch_extend_width,
					 int N_zones_per_right_angle)
		const;

private:
	// verify that grid spacing evenly divides grid sizes
	// in both directions; no-op if ok, error_exit() if not ok
	void verify_grid_spacing_ok(int N_zones_per_right_angle)
		const;
	};

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif
