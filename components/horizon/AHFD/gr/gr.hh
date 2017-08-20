#ifndef AHFD_GR_GR_H
#define AHFD_GR_GR_H

// gr.hh -- header file for general relativity code
// $Header$

// everything in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

//
// number of spatial dimensions in the main Cactus grid
// and in our trial-horizon-surface grid
//
enum	{ N_GRID_DIMS = 3, N_HORIZON_DIMS = 2 };

//
// this enum specifies what kind of surface we want
//
enum    a_surface_definition
        {
        definition_error,      // this value is illegal
        definition_expansion,  // apparent horizon
        definition_inner_expansion, // expansion Theta_(l), ingoing null normal
        definition_mean_curvature, // mean curvature
        definition_expansion_product // product of Theta_(n) and Theta_(l)
        };

//
// this enum specifies how the surface definition is modified
//
enum    a_surface_modification
        {
        modification_error,     // this value is illegal
        modification_none,      // no modification
        modification_radius,    // times coordinate radius
        modification_radius2    // times coordinate radius^2
#if 0
        modification_mean_radius, // times mean coordinate radius
        modification_areal_radius // times areal radius
#endif
        };

//
// this enum specifies how we select the surface
//
enum    a_surface_selection
        {
        selection_error,        // this value is illegal
        selection_definition,   // use the surface's definition
        selection_mean_coordinate_radius, // mean coordinate radius (cheap)
        selection_areal_radius,  // areal radius
        selection_expansion_mean_coordinate_radius, // expansion times mean coordinate radius
        selection_expansion_areal_radius // expansion times areal radius
        };

//
// this struct specifies what to calculate
//
struct  what_to_compute
        {
        // how Theta is calculated
        a_surface_definition surface_definition;
        // how Theta is modified
        a_surface_modification surface_modification;
        // what is solved for
        a_surface_selection surface_selection;
        // the desired value (expansion, areal radius, etc.)
        fp desired_value;
        
        what_to_compute ()
          : surface_definition (definition_error),
            surface_modification (modification_error),
            surface_selection (selection_error),
            desired_value (0.0)
        { }
        };

//
// this enum holds the (a) decoded  Jacobian_compute_method  parameter,
// i.e. it specifies how we compute the (a) Jacobian matrix
//
enum	Jacobian_compute_method
	{
	Jacobian__numerical_perturbation,
	Jacobian__symbolic_diff_with_FD_dr,
	Jacobian__symbolic_diff // no comma
	};

//
// this enum describes the status of an  expansion()  or  expansion_Jacobian()
// computation which returns (i.e. which doesn't abort the Cactus run)
//
enum	expansion_status
	{
	// successful computation
	// ... this is also returned for dummy computations
	expansion_success,

	// non-finite() (i.e. +/-infinity or NaN) values found in h
	expansion_failure__surface_nonfinite,

	// surface is too large
	// ... this value is never returned by  expansion()  or
	//      expansion_Jacobian() , but is placed in this enum
	//     for the convenience of higher-level software which
	//     wishes to have such an error condition
	expansion_failure__surface_too_large,

	// geometry interpolator returns CCTK_ERROR_INTERP_POINT_OUTSIDE
	expansion_failure__surface_outside_grid,

	// geometry interpolator returns CCTK_ERROR_INTERP_POINT_EXCISED
	expansion_failure__surface_in_excised_region,	// (not implemented yet)

	// non-finite (i.e. +/-infinity or NaN) values found
	// in interpolated geometry (g_ij, K_ij, partial_k g_ij)
	expansion_failure__geometry_nonfinite,

	// interpolated g_ij isn't positive definite
	expansion_failure__gij_not_positive_definite // no comma
	};

//******************************************************************************

//
// This structure holds all the information we need about the Cactus grid
// and gridfns in order to interpolate the geometry information.
//
struct	cactus_grid_info
	{
	int coord_system_handle;	// Cactus coordinate system handle

	// this describes the semantics of the Cactus gij gridfns:
	// true ==> the Cactus g_ij are conformal values as per
	//          CactusEinstein/StaticConformal and the  psi  gridfn
	// false ==> the Cactus g_ij are the physical metric
	bool use_Cactus_conformal_metric;

	// Cactus variable indices of geometry variables
        int mask_varindex;
	int g_dd_11_varindex, g_dd_12_varindex, g_dd_13_varindex,
			      g_dd_22_varindex, g_dd_23_varindex,
						g_dd_33_varindex;
	int K_dd_11_varindex, K_dd_12_varindex, K_dd_13_varindex,
			      K_dd_22_varindex, K_dd_23_varindex,
						K_dd_33_varindex;
	int psi_varindex;	// unused if !use_Cactus_conformal_metric
	};

//
// This structure holds information for computing the spacetime geometry.
// This is normally done by interpolating $g_{ij}$ and $K_{ij}$ from the
// usual Cactus grid, but can optionally instead by done by hard-wiring
// the Schwarzschild geometry in Eddington-Finkelstein coordinates.
//
struct	geometry_info
	{
	// normally false; true for testing/debugging
	bool hardwire_Schwarzschild_EF_geometry;

	// parameters for the normal case
	//   hardwire_Schwarzschild_EF_geometry == false
	int operator_handle;		// Cactus handle to interpolation op
	int param_table_handle;		// Cactus handle to parameter table
					// for the interpolator

	// parameters for  hardwire_Schwarzschild_EF_geometry == true
	fp geometry__Schwarzschild_EF__x_posn;	// x posn of Schwarzschild BH
	fp geometry__Schwarzschild_EF__y_posn;	// y posn of Schwarzschild BH
	fp geometry__Schwarzschild_EF__z_posn;	// z posn of Schwarzschild BH
	fp geometry__Schwarzschild_EF__mass;	// mass of Schwarzschild BH
	fp geometry__Schwarzschild_EF__epsilon;	// use z axis limits if
						// (x^2+y^2)/r^2 <= this
	fp geometry__Schwarzschild_EF__Delta_xyz;// "grid spacing" for FD
						// computation of partial_k g_ij

	// should we check various gridfns for NaNs?
	bool check_that_h_is_finite;
	bool check_that_geometry_is_finite;

	// horizon may shrink?
	bool mask_is_noshrink;
	};

//
// This struct holds parameters for computing the Jacobian matrix.
//
struct	Jacobian_info
	{
	enum Jacobian_compute_method     Jacobian_compute_method;
	enum Jacobian_store_solve_method Jacobian_store_solve_method;
	fp perturbation_amplitude;
	};

//
// This struct holds information on what to do if various error/warning
// conditions occur.
//
struct	error_info
	{
	// CCTK_VWarn() level for point outside (or too close to any
	// boundary of) the Cactus grid, i.e. geometry interpolator returns
	// CCTK_ERROR_INTERP_POINT_OUTSIDE
	// ... warning level if error occurs on first Newton iteration,
	//     i.e. when evaluating expansion/Jacobian for initial guess
	int warn_level__point_outside__initial;
	// ... warning level if error occurs on a subsequent Newton iteration
	int warn_level__point_outside__subsequent;

	// CCTK_VWarn() level for the user set
	//  check_that_geometry_is_finite = "true"
	// but the Cactus configure process failed to find the finite()
	// function ==> we have to skip the check
	int warn_level__skipping_finite_check;

	// CCTK_VWarn() level for infinity or NaN in interpolated geometry
	// (g_ij, partial_k g_ij, and/or K_ij)
	int warn_level__nonfinite_geometry;

	// CCTK_VWarn() level for interpolated g_ij isn't positive definite
	// (usually this means we're too close to a singularity)
	// ... warning level if error occurs on first Newton iteration,
	//     i.e. when evaluating expansion/Jacobian for initial guess
	int warn_level__gij_not_positive_definite__initial;
	// ... warning level if error occurs on a subsequent Newton iteration
	int warn_level__gij_not_positive_definite__subsequent;
	};

//******************************************************************************

//
// prototypes for functions visible outside their source files
//

// expansion.cc
enum expansion_status
  expansion(patch_system* ps_ptr,
            const struct what_to_compute& comput_info,
	    const struct cactus_grid_info& cgi,
	    const struct geometry_info& gi,
	    const struct error_info& error_info, bool initial_flag,
	    bool Jacobian_flag = false,
	    bool print_msg_flag = false,
	    jtutil::norm<fp>* H_norms_ptr = NULL,
	    jtutil::norm<fp>* expansion_H_norms_ptr = NULL,
	    jtutil::norm<fp>* inner_expansion_H_norms_ptr = NULL,
	    jtutil::norm<fp>* product_expansion_H_norms_ptr = NULL,
	    jtutil::norm<fp>* mean_curvature_H_norms_ptr = NULL);

// expansion_Jacobian.cc
enum expansion_status
  expansion_Jacobian(patch_system* ps_ptr, Jacobian* Jac_ptr,
                     const struct what_to_compute& comput_info,
		     const struct cactus_grid_info& cgi,
		     const struct geometry_info& gi,
		     const struct Jacobian_info& Jacobian_info,
		     const struct error_info& error_info, bool initial_flag,
		     bool print_msg_flag = false);

// Schwarzschild_EF.cc
void Schwarzschild_EF_geometry(patch_system& ps,
			       const struct geometry_info& gi,
			       bool msg_flag);

// misc-gr.cc
enum Jacobian_compute_method
  decode_Jacobian_compute_method(const char Jacobian_compute_method_string[]);
const char* expansion_status_string(enum expansion_status status);

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif
