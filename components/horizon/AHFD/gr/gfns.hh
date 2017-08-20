#ifndef AHFD_GR_GFNS_H
#define AHFD_GR_GFNS_H

// gfns.hh -- define gfns of all gridfns
// $Header$

// everything in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

//
// For a skeletal patch system we only use the ghosted gridfns.
//

namespace gfns
	  {

// ghosted gridfns
enum	{
	ghosted_min_gfn = -3,		// must set this by hand so
					// ghosted_max_gfn is still < 0
	gfn__h = ghosted_min_gfn,
	gfn__save_h,
	ghosted_max_gfn = gfn__save_h
	};

// nominal gridfns
enum	{
	nominal_min_gfn = 1,

	//
	// for a skeletal patch system we don't need any nominal gridfns
	//
	skeletal_nominal_max_gfn = nominal_min_gfn - 1,

	//
	// most of these gridfns have access macros in "cg.hh";
	// the ones that don't are marked explicitly
	//
	gfn__global_x = nominal_min_gfn,	// no access macro
	gfn__global_y,				// no access macro
	gfn__global_z,				// no access macro
        
	gfn__global_xx,				// no access macro
	gfn__global_xy,				// no access macro
	gfn__global_xz,				// no access macro
	gfn__global_yy,				// no access macro
	gfn__global_yz,				// no access macro
	gfn__global_zz,				// no access macro

        gfn__mask,				// no access macro
        gfn__partial_d_mask_1,			// no access macro
        gfn__partial_d_mask_2,			// no access macro
        gfn__partial_d_mask_3,			// no access macro

	gfn__g_dd_11,
	gfn__g_dd_12,
	gfn__g_dd_13,
	gfn__g_dd_22,
	gfn__g_dd_23,
	gfn__g_dd_33,
	gfn__partial_d_g_dd_111,
	gfn__partial_d_g_dd_112,
	gfn__partial_d_g_dd_113,
	gfn__partial_d_g_dd_122,
	gfn__partial_d_g_dd_123,
	gfn__partial_d_g_dd_133,
	gfn__partial_d_g_dd_211,
	gfn__partial_d_g_dd_212,
	gfn__partial_d_g_dd_213,
	gfn__partial_d_g_dd_222,
	gfn__partial_d_g_dd_223,
	gfn__partial_d_g_dd_233,
	gfn__partial_d_g_dd_311,
	gfn__partial_d_g_dd_312,
	gfn__partial_d_g_dd_313,
	gfn__partial_d_g_dd_322,
	gfn__partial_d_g_dd_323,
	gfn__partial_d_g_dd_333,
	gfn__K_dd_11,
	gfn__K_dd_12,
	gfn__K_dd_13,
	gfn__K_dd_22,
	gfn__K_dd_23,
	gfn__K_dd_33,

	gfn__psi,				// no access macro
	gfn__partial_d_psi_1,			// no access macro
	gfn__partial_d_psi_2,			// no access macro
	gfn__partial_d_psi_3,			// no access macro

	gfn__mean_curvature,			// no access macro
	gfn__save_mean_curvature,		// no access macro

	gfn__Theta,
	gfn__partial_Theta_wrt_partial_d_h_1,
	gfn__partial_Theta_wrt_partial_d_h_2,
	gfn__partial_Theta_wrt_partial_dd_h_11,
	gfn__partial_Theta_wrt_partial_dd_h_12,
	gfn__partial_Theta_wrt_partial_dd_h_22,
	gfn__Delta_h,
	gfn__save_Theta,
	gfn__old_Theta,
	gfn__zero,
	gfn__one,
	nominal_max_gfn = gfn__one		// no comma
	};

	  }	// namespace gfns::

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif
