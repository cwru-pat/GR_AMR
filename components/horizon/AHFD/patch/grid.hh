// grid.hh -- classes for a 2D uniform tensor-product grid
// $Header$
//
#ifndef AHFD_PATCH_GRID_H
#define AHFD_PATCH_GRID_H

#include <cstdio>
#include <assert.h>
#include <math.h>
#include "../cctk.h"
#include "../jtutil/util.hh"
#include "../jtutil/array.hh"
#include "../jtutil/linear_map.hh"
#include "coords.hh"
// grid_arrays - data arrays for a 2D tensor-product grid
// grid - uniform 2D tensor-product grid
//

//
// prerequisites:
//    <stdio.h>
//    <assert.h>
//    <math.h>
//    "cctk.h" or "fake_cctk.h"
//    "stdc.h"
//    "config.hh"
//    "../jtutil/util.hh"
//    "../jtutil/array.hh"
//    "../jtutil/linear_map.hh"
//    "coords.hh"
//


// everything in this file is inside this namespace
namespace AHFinderDirect
	  {

//*****************************************************************************

//
// grid_arrays - data arrays for a 2D tensor-product grid
//
// This is a helper class for class grid (below).  This class stores
// most of the actual grid function (gridfn) data arrays for a uniform
// tensor-product 2D grid.
//
// The integer grid coordinates are (irho,isigma).  This class deals
// with the grid solely at the level of arrays with integer subscripts;
// the derived class  grid  deals with the floating-point coordinates
// related to those subscripts.
//
// The grid has a nominal extent, surrounded by "ghost zones" on each
// side for finite differencing purposes.
//
// There are separate sets of nominal-grid and ghosted-grid gridfns.
// We identify a gridfn by a small-integer "grid function number", a.k.a.
// "gfn".  There are separate gfns for nominal and ghosted gridfns.
// In a very few places we refer to "unknown-grid" gridfns; these might
// be either nominal-grid or ghosted-grid.
//
// For our application (apparent horizon finding), it's useful for the
// storage for a single gridfn to be contiguous *across all patches*.
// (Note this means that the set of all our gridfns is *not* contiguous!)
// To accomplish this, we don't allocate the gridfns when we're created,
// but rather later, with a separate call  setup_gridfn_storage() .
// This way higher-level code can first create all patches, then count
// the total amount of storage used, allocate it, then finally call each
// patch again to set up its gridfns appropriately.
//

class	grid_arrays
	{
public:
	//
	// ***** {min,max}_{rho,sigma} "sides" of grid *****
	//

	//
	// A grid has 4 (angular) "sides", which we identify as
	// {min,max}_{rho,sigma}.  Given a side, we define coordinates
	// (perpendicular,parallel) to it, normally abbreviated to
	// (perp,par).
	//
	// As well as functions directly referring to a specific side,
	// we also support referring to one of these chosen at run-time,
	// via Boolean flags:
	//
	//	// generic (irho,isigma) coordinate
	//	iang = want_rho ? irho : isigma
	//
	//	// opposite (irho,isigma) coordinate
	//	ixang = want_rho ? isigma : irho
	//
	//	// generic (min,max) direction
	//	minmax = want_min ? min : max
	//
	// FIXME: This system of Boolean flags works ok, but it requires
	//	  a lot of repetitive code conditional-expression functions
	//	  in this class.  Is there a cleaner solution?

	// there are precisely this many possible sides
	enum { N_sides = 4 };

	// we specify {min,max} with a Boolean  want_min
	// ... values for want_min
	//     FIXME: these should really be bool, but then we couldn't
	//            use the "enum hack" for in-class constants
	enum { side_is_min = true, side_is_max = false };

	// we specify {rho,sigma} with a Boolean  want_rho
	// ... values for wanr_rho
	//     FIXME: these should really be bool, but then we couldn't
	//            use the "enum hack" for in-class constants
	enum { side_is_rho = true, side_is_sigma = false };

	// human-readable names for the sides (for debugging)
	static const char* minmax_name(bool minmax)
		{ return minmax ? "min" : "max"; }
	static const char* iang_name(bool want_rho)
		{ return want_rho ? "irho" : "isigma"; }


	//
	// ***** array info *****
	//
public:

	// nominal-grid min/max/sizes
	int min_irho() const { return min_irho_; }
	int max_irho() const { return max_irho_; }
	int min_isigma() const { return min_isigma_; }
	int max_isigma() const { return max_isigma_; }
	int min_iang(bool want_rho) const
		{ return want_rho ? min_irho() : min_isigma(); }
	int max_iang(bool want_rho) const
		{ return want_rho ? max_irho() : max_isigma(); }
	int minmax_iang(bool want_min, bool want_rho) const
		{ return want_min ? min_iang(want_rho) : max_iang(want_rho); }
	int N_irho() const
		{ return jtutil::how_many_in_range(min_irho(), max_irho()); }
	int N_isigma() const
		{
		return jtutil::how_many_in_range(min_isigma(), max_isigma());
		}
	int N_grid_points() const
		{ return N_irho() * N_isigma(); }

	// ghosted-grid min/max/sizes
	int ghosted_min_irho() const { return ghosted_min_irho_; }
	int ghosted_max_irho() const { return ghosted_max_irho_; }
	int ghosted_min_isigma() const
		{ return ghosted_min_isigma_; }
	int ghosted_max_isigma() const
		{ return ghosted_max_isigma_; }
	int ghosted_min_iang(bool want_rho) const
		{
		return want_rho ? ghosted_min_irho()
				: ghosted_min_isigma();
		}
	int ghosted_max_iang(bool want_rho) const
		{
		return want_rho ? ghosted_max_irho()
				: ghosted_max_isigma();
		}
	int ghosted_minmax_iang(bool want_min, bool want_rho) const
		{
		return want_min ? ghosted_min_iang(want_rho)
				: ghosted_max_iang(want_rho);
		}
	int ghosted_N_irho() const
		{
		return jtutil::how_many_in_range(ghosted_min_irho(),
						 ghosted_max_irho());
		}
	int ghosted_N_isigma() const
		{
		return jtutil::how_many_in_range(ghosted_min_isigma(),
						 ghosted_max_isigma());
		}
	int ghosted_N_grid_points() const
		{ return ghosted_N_irho() * ghosted_N_isigma(); }

	// "effective" grid min/max/sizes
	// (= dynamic select between nominal and full grids)
	int effective_min_irho(bool want_ghost_zones) const
		{
		return want_ghost_zones ? ghosted_min_irho() : min_irho();
		}
	int effective_max_irho(bool want_ghost_zones) const
		{
		return want_ghost_zones ? ghosted_max_irho() : max_irho();
		}
	int effective_min_isigma(bool want_ghost_zones) const
		{
		return want_ghost_zones ? ghosted_min_isigma() : min_isigma();
		}
	int effective_max_isigma(bool want_ghost_zones) const
		{
		return want_ghost_zones ? ghosted_max_isigma() : max_isigma();
		}
	int effective_N_irho(bool want_ghost_zones) const
		{ return want_ghost_zones ? ghosted_N_irho() : N_irho(); }
	int effective_N_isigma(bool want_ghost_zones) const
		{ return want_ghost_zones ? ghosted_N_isigma() : N_isigma(); }


	//
	// ***** ghost zones *****
	//
public:

	// ghost zone min/max perpendicular coordinates
	int min_rho_ghost_zone__min_iperp() const
		{ return ghosted_min_irho(); }
	int min_rho_ghost_zone__max_iperp() const
		{ return min_irho() - 1; }
	int max_rho_ghost_zone__min_iperp() const
		{ return max_irho() + 1; }
	int max_rho_ghost_zone__max_iperp() const
		{ return ghosted_max_irho(); }
	int min_sigma_ghost_zone__min_iperp() const
		{ return ghosted_min_isigma(); }
	int min_sigma_ghost_zone__max_iperp() const
		{ return min_isigma() - 1; }
	int max_sigma_ghost_zone__min_iperp() const
		{ return max_isigma() + 1; }
	int max_sigma_ghost_zone__max_iperp() const
		{ return ghosted_max_isigma(); }
	int minmax_ang_ghost_zone__min_iperp(bool want_min, bool want_rho) const
		{
		return want_min
		       ? (want_rho ? min_rho_ghost_zone__min_iperp()
				   : min_sigma_ghost_zone__min_iperp())
		       : (want_rho ? max_rho_ghost_zone__min_iperp()
				   : max_sigma_ghost_zone__min_iperp());
		}
	int minmax_ang_ghost_zone__max_iperp(bool want_min, bool want_rho) const
		{
		return want_min
		       ? (want_rho ? min_rho_ghost_zone__max_iperp()
				   : min_sigma_ghost_zone__max_iperp())
		       : (want_rho ? max_rho_ghost_zone__max_iperp()
				   : max_sigma_ghost_zone__max_iperp());
		}

	// ghost zone min/max parallel coordinates
	// ... not including corners
	int rho_ghost_zone_without_corners__min_ipar() const
		{ return min_isigma(); }
	int rho_ghost_zone_without_corners__max_ipar() const
		{ return max_isigma(); }
	int sigma_ghost_zone_without_corners__min_ipar() const
		{ return min_irho(); }
	int sigma_ghost_zone_without_corners__max_ipar() const
		{ return max_irho(); }
	int ang_ghost_zone_without_corners__min_ipar(bool want_rho) const
		{
		return want_rho ?   rho_ghost_zone_without_corners__min_ipar()
				: sigma_ghost_zone_without_corners__min_ipar();
		}
	int ang_ghost_zone_without_corners__max_ipar(bool want_rho) const
		{
		return want_rho ?   rho_ghost_zone_without_corners__max_ipar()
				: sigma_ghost_zone_without_corners__max_ipar();
		}
	// ... including corners
	int rho_ghost_zone_with_corners__min_ipar() const
		{ return ghosted_min_isigma(); }
	int rho_ghost_zone_with_corners__max_ipar() const
		{ return ghosted_max_isigma(); }
	int sigma_ghost_zone_with_corners__min_ipar() const
		{ return ghosted_min_irho(); }
	int sigma_ghost_zone_with_corners__max_ipar() const
		{ return ghosted_max_irho(); }
	int ang_ghost_zone_with_corners__min_ipar(bool want_rho) const
		{
		return want_rho ?   rho_ghost_zone_with_corners__min_ipar()
				: sigma_ghost_zone_with_corners__min_ipar();
		}
	int ang_ghost_zone_with_corners__max_ipar(bool want_rho) const
		{
		return want_rho ?   rho_ghost_zone_with_corners__max_ipar()
				: sigma_ghost_zone_with_corners__max_ipar();
		}


	//
	// ***** grid-point validity and membership predicates *****
	//
public:
	bool is_valid_irho(int irho) const
		{ return (irho >= min_irho()) && (irho <= max_irho()); }
	bool is_valid_isigma(int isigma) const
		{ return (isigma >= min_isigma()) && (isigma <= max_isigma()); }
	bool is_in_nominal_grid(int irho, int isigma) const
		{ return is_valid_irho(irho) && is_valid_isigma(isigma); }

	bool is_valid_ghosted_irho(int irho) const
		{
		return    (irho >= ghosted_min_irho())
		       && (irho <= ghosted_max_irho());
		}
	bool is_valid_ghosted_isigma(int isigma) const
		{
		return    (isigma >= ghosted_min_isigma())
		       && (isigma <= ghosted_max_isigma());
		}
	bool is_in_ghosted_grid(int irho, int isigma) const
		{
		return    is_valid_ghosted_irho(irho)
		       && is_valid_ghosted_isigma(isigma);
		}

	bool is_in_ghost_zone(int irho, int isigma) const
		{
		return     is_in_ghosted_grid(irho, isigma)
		       && !is_in_nominal_grid(irho, isigma);
		}


	//
	// ***** gfn ranges and validity predicates *****
	//
public:
	// gfn ranges
	int min_gfn() const
		{
		assert(gridfn_data_ != NULL);
		return (*gridfn_data_).min_i();
		}
	int max_gfn() const
		{
		assert(gridfn_data_ != NULL);
		return (*gridfn_data_).max_i();
		}
	int N_gridfns() const
		{ return jtutil::how_many_in_range(min_gfn(), max_gfn()); }
	int ghosted_min_gfn() const
		{
		assert(ghosted_gridfn_data_ != NULL);
		return (*ghosted_gridfn_data_).min_i();
		}
	int ghosted_max_gfn() const
		{
		assert(ghosted_gridfn_data_ != NULL);
		return (*ghosted_gridfn_data_).max_i();
		}
	int ghosted_N_gridfns() const
		{
		return jtutil::how_many_in_range(ghosted_min_gfn(),
						 ghosted_max_gfn());
		}

	// gfn validity predicates
	bool is_valid_gfn(int gfn) const
		{ return (gfn >= min_gfn()) && (gfn <= max_gfn()); }
	bool is_valid_ghosted_gfn(int gfn) const
		{
		return (gfn >= ghosted_min_gfn()) && (gfn <= ghosted_max_gfn());
		}


	//
	// ***** gridfns *****
	//
	// n.b. access to rvalue gridfn data must be via references
	//	in order to allow using  gridfn(...)  as the operand
	//	of a unary & (address-of) operator
	//
public:
	// access to nominal-grid gridfn data
	// ... rvalue
	const fp& gridfn(int gfn,   int irho, int isigma) const
		{
		assert(gridfn_data_ != NULL);
		return (*gridfn_data_)(gfn, irho, isigma);
		}
	// ... lvalue
	      fp& gridfn(int gfn,   int irho, int isigma)
		{
		assert(gridfn_data_ != NULL);
		return (*gridfn_data_)(gfn, irho, isigma);
		}

	// access to ghosted-grid gridfn data
	// ... rvalue
	const fp& ghosted_gridfn(int gfn,   int irho, int isigma) const
		{
		assert(gridfn_data_ != NULL);
		return (*ghosted_gridfn_data_)(gfn, irho, isigma);
		}
	// ... lvalue
	      fp& ghosted_gridfn(int gfn,   int irho, int isigma)
		{
		assert(gridfn_data_ != NULL);
		return (*ghosted_gridfn_data_)(gfn, irho, isigma);
		}

	// access to unknown-grid gridfn data
	// (either nominal or ghosted, depending on Boolean flag)
	// ... rvalue
	const fp& unknown_gridfn(bool ghosted_flag,
				 int unknown_gfn,   int irho, int isigma)
		const
		{
		return ghosted_flag ? ghosted_gridfn(unknown_gfn, irho,isigma)
				    :         gridfn(unknown_gfn, irho,isigma);
		}
	// ... lvalue
	      fp& unknown_gridfn(bool ghosted_flag,
				 int unknown_gfn,   int irho, int isigma)
		{
		return ghosted_flag ? ghosted_gridfn(unknown_gfn, irho,isigma)
				    :         gridfn(unknown_gfn, irho,isigma);
		}

	// subscripting info
	int gfn_stride() const
		{
		assert(gridfn_data_ != NULL);
		return gridfn_data_->subscript_stride_i();
		}
	int irho_stride() const
		{
		assert(gridfn_data_ != NULL);
		return gridfn_data_->subscript_stride_j();
		}
	int isigma_stride() const
		{
		assert(gridfn_data_ != NULL);
		return gridfn_data_->subscript_stride_k();
		}
	int iang_stride(bool want_rho) const
		{ return want_rho ? irho_stride() : isigma_stride(); }
	int ghosted_gfn_stride() const
		{
		assert(ghosted_gridfn_data_ != NULL);
		return ghosted_gridfn_data_->subscript_stride_i();
		}
	int ghosted_irho_stride() const
		{
		assert(ghosted_gridfn_data_ != NULL);
		return ghosted_gridfn_data_->subscript_stride_j();
		}
	int ghosted_isigma_stride() const
		{
		assert(ghosted_gridfn_data_ != NULL);
		return ghosted_gridfn_data_->subscript_stride_k();
		}
	int ghosted_iang_stride(bool want_rho) const
		{
		return want_rho ? ghosted_irho_stride()
				: ghosted_isigma_stride();
		}

	// validity predicates for 1-D 0-origin grid point number (gpn)
	bool is_valid_gpn(int gpn) const
		{ return (gpn >= 0) && (gpn < N_grid_points()); }
	bool is_valid_ghosted_gpn(int gpn) const
		{ return (gpn >= 0) && (gpn < ghosted_N_grid_points()); }

	// convert (irho,isigma) <--> 1-D 0-origin grid point number (gpn)
	int gpn_of_irho_isigma(int irho, int isigma) const
		{
		assert( is_valid_irho(irho) );
		assert( is_valid_isigma(isigma) );
		return   (irho   - min_irho()  ) * irho_stride()
		       + (isigma - min_isigma()) * isigma_stride();
		}
	int ghosted_gpn_of_irho_isigma(int irho, int isigma) const
		{
		assert( is_valid_ghosted_irho(irho) );
		assert( is_valid_ghosted_isigma(isigma) );
		return
		     (irho   - ghosted_min_irho()  ) * ghosted_irho_stride()
		   + (isigma - ghosted_min_isigma()) * ghosted_isigma_stride();
		}
	// ... current implementation assumes (& verifies) isigma is contiguous
	void irho_isigma_of_gpn(int gpn, int& irho, int& isigma) const
		{
		assert( is_valid_gpn(gpn) );
		assert( isigma_stride() == 1 );	// implementation restriction
		irho   = min_irho()   + gpn / N_isigma();
		isigma = min_isigma() + gpn % N_isigma();
		assert( is_valid_irho(irho) );
		assert( is_valid_isigma(isigma) );
		}
	// ... current implementation assumes (& verifies) isigma is contiguous
	void ghosted_irho_isigma_of_gpn(int gpn, int& irho, int& isigma) const
		{
		assert( is_valid_ghosted_gpn(gpn) );
		assert( ghosted_isigma_stride() == 1 );	// implementation
							// restriction
		irho   = ghosted_min_irho()   + gpn / ghosted_N_isigma();
		isigma = ghosted_min_isigma() + gpn % ghosted_N_isigma();
		assert( is_valid_ghosted_irho(irho) );
		assert( is_valid_ghosted_isigma(isigma) );
		}

	// low-level access to data arrays (!!dangerous!!)
	const fp* gridfn_data_array(int gfn) const
		{ return & gridfn(gfn, min_irho(),min_isigma()); }
	      fp* gridfn_data_array(int gfn)
		{ return & gridfn(gfn, min_irho(),min_isigma()); }
	const fp* ghosted_gridfn_data_array(int ghosted_gfn) const
		{
		return & ghosted_gridfn(ghosted_gfn, ghosted_min_irho(),
						     ghosted_min_isigma());
		}
	      fp* ghosted_gridfn_data_array(int ghosted_gfn)
		{
		return & ghosted_gridfn(ghosted_gfn, ghosted_min_irho(),
						     ghosted_min_isigma());
		}


	//
	// ***** argument structures for constructor et al *****
	//
public:
	// these structures bundle related arguments together so we don't
	// have 20+ (!) separate arguments to our top-level constructors
	struct	grid_array_pars
		{
		int min_irho,   max_irho;
		int min_isigma, max_isigma;
		int min_rho_ghost_zone_width,   max_rho_ghost_zone_width;
		int min_sigma_ghost_zone_width, max_sigma_ghost_zone_width;
		};
	struct	gridfn_pars
		{
		int min_gfn, max_gfn;

		// gridfn storage will be automatically allocated
		// if pointer is NULL; any 0 strides are automatically
		// set to C-style row-major subscripting
		fp* storage_array;
		int gfn_stride, irho_stride, isigma_stride;
		};

	//
	// ***** constructor, gridfn setup, destructor *****
	//
public:
	// construct with no gridfns
	grid_arrays(const grid_array_pars& grid_array_pars_in);

	// set up storage for gridfns
	void setup_gridfn_storage(const gridfn_pars& gridfn_pars_in,
				  const gridfn_pars& ghosted_gridfn_pars_in);

	~grid_arrays();

private:
	// we forbid copying and passing by value
	// by declaring the copy constructor and assignment operator
	// private, but never defining them
	grid_arrays(const grid_arrays& rhs);
	grid_arrays& operator=(const grid_arrays& rhs);

private:
	//
	// ***** the actual gridfn storage arrays *****
	//
	// n.b. these pointers are *first* data member in this class
	// ==> possibly slightly faster access (0 offset from pointer)
	// ... indices are (gfn, irho, isigma)
	jtutil::array3d<fp>* gridfn_data_;
	jtutil::array3d<fp>* ghosted_gridfn_data_;

	// gfn bounds
	const int min_gfn_,         max_gfn_;
	const int ghosted_min_gfn_, ghosted_max_gfn_;

	// nominal grid min/max bounds
	const int min_irho_,   max_irho_;
	const int min_isigma_, max_isigma_;

	// full grid min/max bounds
	const int ghosted_min_irho_,   ghosted_max_irho_;
	const int ghosted_min_isigma_, ghosted_max_isigma_;
	};

//******************************************************************************

//
// grid - uniform 2D tensor-product grid
//
// The grid is uniform in the floating point grid coordinates (rho,sigma).
// There is also some (limited) support for expressing these coordinates
// in degrees (drho,dsigma); this is useful for humans trying to specify
// things in parameter files.
//
// The nominal (not including the ghost zones) angular grid boundaries
// may coincide with grid points, or they may be at "half-integer" grid
// coordinates.  That is, suppose we have a unit grid spacing, and a boundary
// at an angular coordinate of 0; then the grid may be either 0, 1, 2, ...,
// or 0.5, 1.5, 2.5, ... .
//

class	grid
	: public grid_arrays
	{
	//
	// ***** low-level access to coordinate maps *****
	//
public:
	// direct (read-only) access to the underlying linear_map objects
	// ... useful for (eg) passing to interpolators
	const jtutil::linear_map<fp>& rho_map() const { return rho_map_; }
	const jtutil::linear_map<fp>& sigma_map() const { return sigma_map_; }
	const jtutil::linear_map<fp>& ang_map(bool want_rho) const
		{ return want_rho ? rho_map() : sigma_map(); }


	//
	// ***** single-axis coordinate conversions *****
	//
public:
	// ... angles in radians
	fp rho_of_irho(int irho) const { return rho_map().fp_of_int(irho); }
	fp sigma_of_isigma(int isigma) const
		{ return sigma_map().fp_of_int(isigma); }
	fp ang_of_iang(bool want_rho, int iang) const
		{
		return want_rho ? rho_of_irho(iang)
				: sigma_of_isigma(iang);
		}

	fp fp_irho_of_rho(fp rho) const
		{ return rho_map().fp_int_of_fp(rho); }
	int irho_of_rho(fp rho, jtutil::linear_map<fp>::noninteger_action
				nia = jtutil::linear_map<fp>::nia_error)
		const
		{ return rho_map().int_of_fp(rho, nia); }
	fp fp_isigma_of_sigma(fp sigma) const
		{ return sigma_map().fp_int_of_fp(sigma); }
	int isigma_of_sigma(fp sigma, jtutil::linear_map<fp>::noninteger_action
				      nia = jtutil::linear_map<fp>::nia_error)
		const
		{ return sigma_map().int_of_fp(sigma, nia); }
	fp fp_iang_of_ang(bool want_rho, fp ang)
		const
		{
		return want_rho ? fp_irho_of_rho(ang)
				: fp_isigma_of_sigma(ang);
		}
	int iang_of_ang(bool want_rho,
			fp ang, jtutil::linear_map<fp>::noninteger_action
				nia = jtutil::linear_map<fp>::nia_error)
		const
		{
		return want_rho ? irho_of_rho(ang, nia)
				: isigma_of_sigma(ang, nia);
		}

	// ... angles in degrees
	fp rho_of_drho(fp drho) const
		{ return jtutil::radians_of_degrees(drho); }
	fp sigma_of_dsigma(fp dsigma) const
		{ return jtutil::radians_of_degrees(dsigma); }
	fp drho_of_rho(fp rho) const
		{ return jtutil::degrees_of_radians(rho); }
	fp dsigma_of_sigma(fp sigma) const
		{ return jtutil::degrees_of_radians(sigma); }
	fp drho_of_irho(int irho) const
		{ return jtutil::degrees_of_radians(rho_of_irho(irho)); }
	fp dsigma_of_isigma(int isigma) const
		{ return jtutil::degrees_of_radians(sigma_of_isigma(isigma)); }

	int irho_of_drho(fp drho, jtutil::linear_map<fp>::noninteger_action
				  nia = jtutil::linear_map<fp>::nia_error)
		const
		{ return irho_of_rho(jtutil::radians_of_degrees(drho), nia); }
	int isigma_of_dsigma(fp dsigma,
			     jtutil::linear_map<fp>::noninteger_action
			        nia = jtutil::linear_map<fp>::nia_error)
		const
		{
		return isigma_of_sigma(jtutil::radians_of_degrees(dsigma), nia);
		}


	//
	// ***** grid info *****
	//
public:
	// grid spacings
	fp delta_rho() const { return rho_map().delta_fp(); }
	fp delta_sigma() const { return sigma_map().delta_fp(); }
	fp delta_drho() const
		{ return jtutil::degrees_of_radians(delta_rho()); }
	fp delta_dsigma() const
		{ return jtutil::degrees_of_radians(delta_sigma()); }
	fp delta_ang(bool want_rho) const
		{ return want_rho ? delta_rho() : delta_sigma(); }
	fp delta_dang(bool want_rho) const
		{ return want_rho ? delta_drho() : delta_dsigma(); }

	// inverse grid spacings
	fp inverse_delta_rho() const { return rho_map().inverse_delta_fp(); }
	fp inverse_delta_sigma() const
		{ return sigma_map().inverse_delta_fp(); }

	// nominal grid min/max
	fp min_rho() const { return min_rho_; }
	fp max_rho() const { return max_rho_; }
	fp min_sigma() const { return min_sigma_; }
	fp max_sigma() const { return max_sigma_; }
	fp minmax_ang(bool want_min, bool want_rho) const
		{
		return want_min ? (want_rho ? min_rho() : min_sigma())
				: (want_rho ? max_rho() : max_sigma());
		}
	fp min_drho() const { return jtutil::degrees_of_radians(min_rho()); }
	fp max_drho() const { return jtutil::degrees_of_radians(max_rho()); }
	fp min_dsigma() const
		{ return jtutil::degrees_of_radians(min_sigma()); }
	fp max_dsigma() const
		{ return jtutil::degrees_of_radians(max_sigma()); }
	fp min_dang(bool want_rho) const
		{ return want_rho ? min_drho() : min_dsigma(); }
	fp max_dang(bool want_rho) const
		{ return want_rho ? max_drho() : max_dsigma(); }

	// ghosted-grid min/max
	fp ghosted_min_rho() const
		{ return rho_of_irho(ghosted_min_irho()); }
	fp ghosted_max_rho() const
		{ return rho_of_irho(ghosted_max_irho()); }
	fp ghosted_min_sigma() const
		{ return sigma_of_isigma(ghosted_min_isigma()); }
	fp ghosted_max_sigma() const
		{ return sigma_of_isigma(ghosted_max_isigma()); }

	// is a given (drho,dsigma) within the grid?
	bool is_valid_drho(fp drho) const
		{
		return    jtutil::fuzzy<fp>::GE(drho, min_drho())
		       && jtutil::fuzzy<fp>::LE(drho, max_drho());
		}
	bool is_valid_dsigma(fp dsigma) const
		{
		return    jtutil::fuzzy<fp>::GE(dsigma, min_dsigma())
		       && jtutil::fuzzy<fp>::LE(dsigma, max_dsigma());
		}

	// reduce a rho/sigma coordinate modulo 2*pi radians (360 degrees)
	// to be within the ghosted grid,
	// or error_exit() if no such value exists
	fp modulo_reduce_rho(fp rho_in) const
		{
		return local_coords
		       ::modulo_reduce_ang(rho_in, ghosted_min_rho(),
						   ghosted_max_rho());
		}
	fp modulo_reduce_sigma(fp sigma_in) const
		{
		return local_coords
		       ::modulo_reduce_ang(sigma_in, ghosted_min_sigma(),
						     ghosted_max_sigma());
		}
	fp modulo_reduce_ang(bool want_rho, fp ang_in) const
		{
		return want_rho ? modulo_reduce_rho(ang_in)
				: modulo_reduce_sigma(ang_in);
		}

	//
	// ***** misc stuff *****
	//
public:
	// human-readable names for the sides (for debugging)
	static const char* ang_name(bool want_rho)
		{ return want_rho ?  "rho" :  "sigma"; }
	static const char* dang_name(bool want_rho)
		{ return want_rho ? "drho" : "dsigma"; }


	//
	// ***** argument structure for constructor *****
	//

	// this structure bundles related arguments together so we don't
	// have 20+ (!) separate arguments to our top-level constructors
	struct	grid_pars		// *** note angles in degrees ***
		{
		fp min_drho,   delta_drho,   max_drho;
		fp min_dsigma, delta_dsigma, max_dsigma;
		};


	//
	// ***** constructor, destructor *****
	//
	grid(const grid_array_pars& grid_array_pars_in,
	     const grid_pars& grid_pars_in);
	// compiler-generated default destructor is ok

private:
	// we forbid copying and passing by value
	// by declaring the copy constructor and assignment operator
	// private, but never defining them
	grid(const grid& rhs);
	grid& operator=(const grid& rhs);

private:
	// range of these is the full grid (including ghost zones)
	const jtutil::linear_map<fp> rho_map_, sigma_map_;

	// angular boundaries of nominal grid
	const fp min_rho_,   max_rho_;
	const fp min_sigma_, max_sigma_;
	};

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif
