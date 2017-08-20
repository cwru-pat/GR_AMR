#ifndef AHFD_PATCH_PATCH_SYSTEM_INFO_H
#define AHFD_PATCH_PATCH_SYSTEM_INFO_H
// patch_system_info.hh -- static info describing various types of patch systems
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
//    "patch_info.hh"
//

// everything in this file is inside this namespace
namespace AHFD
	  {

//******************************************************************************

//
// This namespace contains static data describing the patch sizes and
// shapes for each type of patch system.  Since this data only describes
// the patch sizes/shapes, we don't distinguish between the different
// boundary conditions.
//

namespace patch_system_info
{
//
// full-sphere patch system
// ... covers all 4pi steradians
//
namespace full_sphere
    {
    enum {
	 patch_number__pz = 0,
	 patch_number__px,
	 patch_number__py,
	 patch_number__mx,
	 patch_number__my,
	 patch_number__mz,
	 N_patches // no comma
	 };
    static const struct patch_info patch_info_array[N_patches]
      = {
	// +z patch (90 x 90 degrees): dmu [ -45,    45], dnu  [ -45,    45]
	 {"+z", patch::patch_is_plus,  'z',  -45.0,  45.0,       -45.0,  45.0},

	// +x patch (90 x 90 degrees): dnu [  45,   135], dphi [ -45,    45]
	 {"+x", patch::patch_is_plus,  'x',   45.0, 135.0,       -45.0,  45.0},

	// +y patch (90 x 90 degrees): dmu [  45,   135], dphi [  45,   135]
	 {"+y", patch::patch_is_plus,  'y',   45.0, 135.0,        45.0, 135.0},

	// -x patch (90 x 90 degrees): dnu [-135,   -45], dphi [ 135,   225]
	 {"-x", patch::patch_is_minus, 'x', -135.0, -45.0,       135.0, 225.0},

	// -y patch (90 x 90 degrees): dmu [-135,   -45], dphi [-135,   -45]
	 {"-y", patch::patch_is_minus, 'y', -135.0, -45.0,      -135.0, -45.0},

	// -z patch (90 x 90 degrees): dmu [ 135,   225], dnu  [ 135,   225]
	 {"-z", patch::patch_is_minus, 'z',  135.0, 225.0,       135.0, 225.0},
	};
    }	// namespace patch_system_info::full_sphere

//
// +z hemisphere (half) patch system
// ... mirror symmetry across z=0 plane
//
namespace plus_z_hemisphere
    {
    enum {
	 patch_number__pz = 0,
	 patch_number__px,
	 patch_number__py,
	 patch_number__mx,
	 patch_number__my,
	 N_patches // no comma
	 };
    static const struct patch_info patch_info_array[N_patches]
      = {
	// +z patch (90 x 90 degrees): dmu [ -45,    45], dnu  [ -45,    45]
	 {"+z", patch::patch_is_plus,  'z',  -45.0,  45.0,       -45.0,  45.0},

	// +x patch (45 x 90 degrees): dnu [  45,    90], dphi [ -45,    45]
	 {"+x", patch::patch_is_plus,  'x',   45.0,  90.0,       -45.0,  45.0},

	// +y patch (45 x 90 degrees): dmu [  45,    90], dphi [  45,   135]
	 {"+y", patch::patch_is_plus,  'y',   45.0,  90.0,        45.0, 135.0},

	// -x patch (45 x 90 degrees): dnu [ -90,   -45], dphi [ 135,   225]
	 {"-x", patch::patch_is_minus, 'x',  -90.0, -45.0,       135.0, 225.0},

	// -y patch (45 x 90 degrees): dmu [ -90,   -45], dphi [-135,   -45]
	 {"-y", patch::patch_is_minus, 'y',  -90.0, -45.0,      -135.0, -45.0},
	};
    }	// namespace patch_system_info::plus_z_hemisphere

//
// +[xy] "vertical" quarter-grid (quadrant) patch system
// two types of boundary conditions:
// ... mirror symmetry across x=0 and y=0 planes
// ... 90 degree periodic rotation symmetry about z axis
//
namespace plus_xy_quadrant
    {
    enum {
	 patch_number__pz = 0,
	 patch_number__px,
	 patch_number__py,
	 patch_number__mz,
	 N_patches // no comma
	 };
    static const struct patch_info patch_info_array[N_patches]
      = {
	// +z patch (45 x 45 degrees): dmu [   0,    45], dnu  [   0,    45]
	 {"+z", patch::patch_is_plus,  'z',  0.0,    45.0,         0.0,  45.0},

	// +x patch (90 x 45 degrees): dnu [  45,   135], dphi [   0,    45]
	 {"+x", patch::patch_is_plus,  'x',  45.0, 135.0,         0.0,  45.0},

	// +y patch (90 x 45 degrees): dmu [  45,   135], dphi [  45,    90]
	 {"+y", patch::patch_is_plus,  'y',  45.0, 135.0,        45.0,  90.0},

	// -z patch (45 x 45 degrees): dmu [ 135,   180], dnu  [ 135,   180]
	 {"-z", patch::patch_is_minus,  'z', 135.0, 180.0,       135.0, 180.0},
	};
    }	// namespace patch_system_info::plus_xy_quadrant

//
// +[xz] "horizontal" quarter-grid (quadrant) patch system
// two types of boundary conditions
// ... mirror symmetry across x=0 plane, z=0 plane
// ... 180 degree periodic rotation symmetry about z axis,
//     mirror symmetry across z=0 plane
//
namespace plus_xz_quadrant
    {
    enum {
	 patch_number__pz = 0,
	 patch_number__px,
	 patch_number__py,
	 patch_number__my,
	 N_patches // no comma
	 };
    static const struct patch_info patch_info_array[N_patches]
      = {
	// +z patch (90 x 45 degrees): dmu [ -45,    45], dnu  [   0,    45]
	 {"+z", patch::patch_is_plus,  'z',  -45.0,  45.0,         0.0,  45.0},

	// +x patch (45 x 90 degrees): dnu [  45,    90], dphi [ -45,    45]
	 {"+x", patch::patch_is_plus,  'x',   45.0,  90.0,       -45.0,  45.0},

	// +y patch (45 x 45 degrees): dmu [  45,    90], dphi [  45,    90]
	 {"+y", patch::patch_is_plus,  'y',   45.0,  90.0,        45.0,  90.0},

	// -y patch (45 x 45 degrees): dmu [ -90,   -45], dphi [ -90,   -45]
	 {"-y", patch::patch_is_minus, 'y',  -90.0, -45.0,       -90.0, -45.0},
	};
    }	// namespace patch_system_info::plus_xz_quadrant_rotating

//
// +[xyz] (octant) patch system
// two types of boundary conditions:
// ... mirror symmetry across x=0 plane, y=0 plane, z=0 plane
// ... 90 degree periodic rotation symmetry about z axis,
//     mirror symmetry across z=0 plane
//
namespace plus_xyz_octant
    {
    enum {
	 patch_number__pz = 0,
	 patch_number__px,
	 patch_number__py,
	 N_patches // no comma
	 };
    static const struct patch_info patch_info_array[N_patches]
      = {
	// +z patch (45 x 45 degrees): dmu [   0,    45], dnu  [   0,    45]
	 {"+z", patch::patch_is_plus,  'z',    0.0,  45.0,         0.0,  45.0},

	// +x patch (45 x 45 degrees): dnu [  45,    90], dphi [   0,    45]
	 {"+x", patch::patch_is_plus,  'x',   45.0,  90.0,         0.0,  45.0},

	// +y patch (45 x 45 degrees): dmu [  45,    90], dphi [  45,    90]
	 {"+y", patch::patch_is_plus,  'y',   45.0,  90.0,        45.0,  90.0},
	};
    }	// namespace patch_system_info::octant_mirrored

	  }	// namespace patch_system_info::

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif
