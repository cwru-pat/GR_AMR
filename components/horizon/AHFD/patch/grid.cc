// grid.cc -- classes for a 2D uniform tensor-product grid
// $Header$
//
// grid_arrays::grid_arrays
// grid_arrays::setup_gridfn_storage
// grid_arrays::~grid_arrays
//
// grid::grid
//

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "../cctk.h"

#include "../jtutil/util.hh"
#include "../jtutil/array.hh"
#include "../jtutil/linear_map.hh"

#include "coords.hh"
#include "grid.hh"

// all the code in this file is inside this namespace
namespace AHFinderDirect
	  {

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

//
// This function constructs a  grid_arrays  object.
//
grid_arrays::grid_arrays(const grid_array_pars& grid_array_pars_in)

	: gridfn_data_(NULL),
	  ghosted_gridfn_data_(NULL),

	  // these are all set properly by setup_gridfn_storage()
	  min_gfn_(0), max_gfn_(0),
	  ghosted_min_gfn_(0), ghosted_max_gfn_(0),

	  min_irho_(grid_array_pars_in.min_irho),
	  max_irho_(grid_array_pars_in.max_irho),
	  min_isigma_(grid_array_pars_in.min_isigma),
	  max_isigma_(grid_array_pars_in.max_isigma),

	  ghosted_min_irho_(grid_array_pars_in.min_irho
			    - grid_array_pars_in.min_rho_ghost_zone_width),
	  ghosted_max_irho_(grid_array_pars_in.max_irho
			    + grid_array_pars_in.max_rho_ghost_zone_width),
	  ghosted_min_isigma_(grid_array_pars_in.min_isigma
			      - grid_array_pars_in.min_sigma_ghost_zone_width),
	  ghosted_max_isigma_(grid_array_pars_in.max_isigma
			      + grid_array_pars_in.max_sigma_ghost_zone_width)
								     // no comma
{ }

//*****************************************************************************

//
// This function sets up the gridfn storage arrays in a  grid_arrays  object.
//
void grid_arrays::setup_gridfn_storage
	(const gridfn_pars& gridfn_pars_in,
	 const gridfn_pars& ghosted_gridfn_pars_in)
{
assert(gridfn_data_ == NULL);
gridfn_data_ = new jtutil::array3d<fp>(gridfn_pars_in.min_gfn,
				       gridfn_pars_in.max_gfn,
				       min_irho(), max_irho(),
				       min_isigma(), max_isigma(),
				       gridfn_pars_in.storage_array,
				       gridfn_pars_in.gfn_stride,
				       gridfn_pars_in.irho_stride,
				       gridfn_pars_in.isigma_stride);

assert(ghosted_gridfn_data_ == NULL);
ghosted_gridfn_data_
	= new jtutil::array3d<fp>
			(ghosted_gridfn_pars_in.min_gfn,
			 ghosted_gridfn_pars_in.max_gfn,
			 ghosted_min_irho(), ghosted_max_irho(),
			 ghosted_min_isigma(), ghosted_max_isigma(),
			 ghosted_gridfn_pars_in.storage_array,
			 ghosted_gridfn_pars_in.gfn_stride,
			 ghosted_gridfn_pars_in.irho_stride,
			 ghosted_gridfn_pars_in.isigma_stride);
}

//******************************************************************************

//
// This function destroys a  grid_arrays  object.
//
grid_arrays::~grid_arrays()
{
delete ghosted_gridfn_data_;
delete gridfn_data_;
}

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

//
// This function constructs a  grid  object.
//
grid::grid(const grid_array_pars& grid_array_pars_in,
	   const grid_pars& grid_pars_in)

	: grid_arrays(grid_array_pars_in),

	  rho_map_(grid_array_pars_in.min_irho
		    - grid_array_pars_in.min_rho_ghost_zone_width,
		   grid_array_pars_in.max_irho
		    + grid_array_pars_in.max_rho_ghost_zone_width,
		   jtutil::radians_of_degrees(
		      grid_pars_in.min_drho
		       - grid_array_pars_in.min_rho_ghost_zone_width
			 * grid_pars_in.delta_drho
					     ),
		   jtutil::radians_of_degrees(grid_pars_in.delta_drho),
		   jtutil::radians_of_degrees(
		      grid_pars_in.max_drho
		       + grid_array_pars_in.max_rho_ghost_zone_width
			 * grid_pars_in.delta_drho
					     )),

	  sigma_map_(grid_array_pars_in.min_isigma
		      - grid_array_pars_in.min_sigma_ghost_zone_width,
		     grid_array_pars_in.max_isigma
		      + grid_array_pars_in.max_sigma_ghost_zone_width,
		     jtutil::radians_of_degrees(
			grid_pars_in.min_dsigma
			 - grid_array_pars_in.min_sigma_ghost_zone_width
			   * grid_pars_in.delta_dsigma
					       ),
		     jtutil::radians_of_degrees(grid_pars_in.delta_dsigma),
		     jtutil::radians_of_degrees(
			grid_pars_in.max_dsigma
			 + grid_array_pars_in.max_sigma_ghost_zone_width
			   * grid_pars_in.delta_dsigma
					       )),

	  min_rho_(jtutil::radians_of_degrees(grid_pars_in.min_drho)),
	  max_rho_(jtutil::radians_of_degrees(grid_pars_in.max_drho)),
	  min_sigma_(jtutil::radians_of_degrees(grid_pars_in.min_dsigma)),
	  max_sigma_(jtutil::radians_of_degrees(grid_pars_in.max_dsigma))
								     // no comma
{ }

//******************************************************************************
//******************************************************************************
//******************************************************************************

	  }	// namespace AHFinderDirect
