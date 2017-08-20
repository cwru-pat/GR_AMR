#ifndef AHFD_DRIVER_BHDIAG_H
#define AHFD_DRIVER_BHDIAG_H
// BH_diagnostics.hh -- header file for BH diagnostics
// $Header$
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <vector>

#include "../jtutil/util_Table.h"
#include "../cctk.h"

#include "../jtutil/util.hh"
#include "../jtutil/array.hh"
#include "../jtutil/cpm_map.hh"
#include "../jtutil/linear_map.hh"

#include "../patch/coords.hh"
#include "../patch/grid.hh"
#include "../patch/fd_grid.hh"
#include "../patch/patch.hh"
#include "../patch/patch_edge.hh"
#include "../patch/patch_interp.hh"
#include "../patch/ghost_zone.hh"
#include "../patch/patch_system.hh"

#include "../elliptic/Jacobian.hh"

#include "../gr/gfns.hh"
#include "../gr/gr.hh"

#include "horizon_sequence.hh"
#include "../AHFD_types.h"

//#include "../AHFD.h"
//
// prerequisites:
//	<stdio.h>
//	"cctk_Arguments.h"
//

// everything in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

//
// This struct holds info for computing black hole diagnostics.
//
struct	BH_diagnostics_info
	{
	enum patch::integration_method integral_method;
	};

//******************************************************************************

//
// A  struct BH_diagnostics  holds all of our black hole diagnostics
// for a single apparent horizon.  These diagnostics are only meaningful
// if the apparent horizon has indeed been found.
//
// Note that all the diagnostics are for the full apparent horizon, even
// if we don't actually store all of it due to patch-system symmetries.
//
struct	BH_diagnostics
	{
public:
	fp origin_x, origin_y, origin_z;
	fp centroid_x, centroid_y, centroid_z;
        fp quadrupole_xx, quadrupole_xy, quadrupole_xz;
        fp quadrupole_yy, quadrupole_yz, quadrupole_zz;
	fp min_radius, max_radius, mean_radius;

	// xyz bounding box
	fp min_x, max_x, min_y, max_y, min_z, max_z;

        fp circumference_xy, circumference_xz, circumference_yz;
        fp area;
        fp expansion;
        fp inner_expansion;
        fp product_expansion;
        fp mean_curvature;

        fp area_gradient;
        fp expansion_gradient;
        fp inner_expansion_gradient;
        fp product_expansion_gradient;
        fp mean_curvature_gradient;
        fp mean_curvature_minimum;
        fp mean_curvature_maximum;
        fp mean_curvature_integral;

public:
	// position of diagnostics in buffer and number of diagnostics
	enum	{
		posn__origin_x = 0, posn__origin_y, posn__origin_z,
		posn__centroid_x, posn__centroid_y, posn__centroid_z,
                posn__quadrupole_xx, posn__quadrupole_xy, posn__quadrupole_xz,
                posn__quadrupole_yy, posn__quadrupole_yz, posn__quadrupole_zz,
		posn__min_radius, posn__max_radius, posn__mean_radius,

		posn__min_x, posn__max_x,
		posn__min_y, posn__max_y,
		posn__min_z, posn__max_z,

		posn__circumference_xy, posn__circumference_xz,
					posn__circumference_yz,
		posn__area,
                posn__expansion,
                posn__inner_expansion,
                posn__product_expansion,
                posn__mean_curvature,

		posn__area_gradient,
                posn__expansion_gradient,
                posn__inner_expansion_gradient,
                posn__product_expansion_gradient,
                posn__mean_curvature_gradient,
                posn__mean_curvature_minimum,
                posn__mean_curvature_maximum,
                posn__mean_curvature_integral,

		N_buffer // no comma	// size of buffer
		};

	// copy diagnostics to/from buffer
	void copy_to_buffer  (      CCTK_REAL buffer[N_buffer]) const;
	void copy_from_buffer(const CCTK_REAL buffer[N_buffer]);

public:
	// compute diagnostics (assuming that apparent horizon has been found)
	void compute(const patch_system& ps,
                     const fp area,
                     const fp mean_expansion,
                     const fp mean_inner_expansion,
                     const fp mean_product_expansion,
                     const fp mean_mean_curvature,
                     const fp area_gradient,
                     const fp mean_expansion_gradient,
                     const fp mean_inner_expansion_gradient,
                     const fp mean_product_expansion_gradient,
                     const fp mean_mean_curvature_gradient,
		     const struct BH_diagnostics_info& BH_diagnostics_info);

	// print (CCTK_VInfo()) a line or two summarizing diagnostics
	void print(int N_horizons, int hn)
		const;

	// create/open output file and write header describing output() fields
	// ... stream is flushed after output to help with
	//     looking at diagnostics while Cactus is still running
	FILE* setup_output_file(const struct IO_info& IO_info,
				int N_horizons, int hn)
		const;

	// output a (long) line of all the diagnostics, to a stdio stream
	// ... stream is flushed after output to help with
	//     looking at diagnostics while Cactus is still running
	void output(FILE* fileptr, const struct IO_info& IO_info)
		const;

	// store the diagnostics in the Cactus variables

	// constructor initializes all diagnostics to 0.0
	BH_diagnostics();

	// no destructor needed, compiler-generated no-op is fine

private:
	// helper function: compute surface integral of specified gridfn
	static
	  fp surface_integral(const patch_system& ps,
			      int src_gfn, bool src_gfn_is_even_across_xy_plane,
					   bool src_gfn_is_even_across_xz_plane,
					   bool src_gfn_is_even_across_yz_plane,
			      enum patch::integration_method method);

private:
	// we forbid copying and passing by value
	// by declaring the copy constructor and assignment operator
	// private, but never defining them
	BH_diagnostics(const BH_diagnostics& rhs);
	BH_diagnostics& operator=(const BH_diagnostics& rhs);
	};

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif
