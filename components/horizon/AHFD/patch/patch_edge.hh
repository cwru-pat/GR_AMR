#ifndef AHFD_PATCH_PATCH_EDGE_H
#define AHFD_PATCH_PATCH_EDGE_H
// patch_edge.hh -- perpendicular/parallel geometry of one side of a patch
// $Header$
//

//
// prerequisites:
//	<stdio.h>
//	<assert.h>
//	<math.h>
//	"stdc.h"
//	"config.hh"
//	"../jtutil/util.hh"
//	"../jtutil/array.hh"
//	"../jtutil/linear_map.hh"
//	"coords.hh"
//	"grid.hh"
//	"fd_grid.hh"
//	"patch.hh"
//

// everything in this file is inside this namespace
namespace AHFD
	  {

//*****************************************************************************

//
// patch_edge -- perpendicular/parallel geometry of one side of a patch
//
// A  patch_edge  object is a very light-weight object which represents
// the basic geometry of a min/max rho/sigma side of a patch, i.e. it
// provides which-side-am-I predicates, coordinate conversions between
// (perp,par) and (rho,sigma), etc.  Every patch has (points to) 4  patch_edge
//  objects, one for each of the patch's sides.  See the comments in
// "patch.hh" for a "big picture" discussion of patches, patch edges,
// ghost zones, and patch interpolation regions.
//
// Note that since  patch_edge  has only  const  member functions
// (and members!), a  patch_edge  object is effectively always  const .
// This means there's no harm in always declaring  patch_edge  objects
// to be  const .
//

class	patch_edge
	{
public:
	//
	// ***** meta-info *****
	//

	// meta-info about patch
	patch& my_patch() const { return my_patch_; }

	// meta-info about edge
	bool is_rho() const { return is_rho_; }
	bool is_min() const { return is_min_; }
	bool perp_is_rho() const { return is_rho(); }
	bool par_is_rho() const { return ! is_rho(); }

	// human-readable {min,max}_{rho,sigma} name (for debugging etc)
	const char* name() const
		{
		return is_min()
		       ? (is_rho() ? "min_rho" : "min_sigma")
		       : (is_rho() ? "max_rho" : "max_sigma");
		}

	// are two edges really the same edge?
	bool operator==(const patch_edge& other_edge) const
		{
		return    ( my_patch() == other_edge.my_patch() )
		       && (   is_rho() == other_edge.  is_rho() )
		       && (   is_min() == other_edge.  is_min() );
		}
	bool operator!=(const patch_edge& other_edge) const
		{ return ! operator==(other_edge); }


	//
	// ***** adjacent edges *****
	//

	// get adjacent edges to our min/max par corners
	const patch_edge& min_par_adjacent_edge() const
		{
		return my_patch()
		       .minmax_ang_patch_edge(grid::side_is_min, par_is_rho());
		}
	const patch_edge& max_par_adjacent_edge() const
		{
		return my_patch()
		       .minmax_ang_patch_edge(grid::side_is_max, par_is_rho());
		}
	const patch_edge& minmax_par_adjacent_edge(bool want_min) const
		{
		return want_min ? min_par_adjacent_edge()
				: max_par_adjacent_edge();
		}


	//
	// ***** gridfn subscripting and coordinate maps *****
	//

	// gridfn strides perpendicular/parallel to the edge
	int perp_stride() const
		{ return my_patch().iang_stride(perp_is_rho()); }
	int par_stride() const
		{ return my_patch().iang_stride(par_is_rho()); }
	int ghosted_perp_stride() const
		{ return my_patch().ghosted_iang_stride(perp_is_rho()); }
	int ghosted_par_stride() const
		{ return my_patch().ghosted_iang_stride(par_is_rho()); }

	// coordinate maps perpendicular/parallel to the edge
	// ... range is that of the grid *including* ghost zones
	const jtutil::linear_map<fp>& perp_map() const
		{ return my_patch().ang_map(perp_is_rho()); }
	const jtutil::linear_map<fp>& par_map() const
		{ return my_patch().ang_map(par_is_rho()); }

	// meta-info about perp/par coordinates
	// ... as (mu,nu,phi) tensor indices
	local_coords::coords_set coords_set_perp() const
		{
		return perp_is_rho() ? my_patch().coords_set_rho()
				     : my_patch().coords_set_sigma();
		}
	local_coords::coords_set coords_set_par() const
		{
		return par_is_rho() ? my_patch().coords_set_rho()
				    : my_patch().coords_set_sigma();
		}


	//
	// ***** coordinate conversions *****
	//

	// coordinate conversions based on ghost zone direction
	// ... (iperp,ipar) <--> (perp,par)
	fp perp_of_iperp(int iperp) const
		{ return my_patch().ang_of_iang(perp_is_rho(), iperp); }
	fp par_of_ipar(int ipar) const
		{ return my_patch().ang_of_iang(par_is_rho(), ipar); }
	fp fp_iperp_of_perp(fp perp) const
		{ return my_patch().fp_iang_of_ang(perp_is_rho(), perp); }
	fp fp_ipar_of_par(fp par) const
		{ return my_patch().fp_iang_of_ang(par_is_rho(), par); }
	int iperp_of_perp(fp perp, jtutil::linear_map<fp>::noninteger_action
				   nia = jtutil::linear_map<fp>::nia_error)
		{ return my_patch().iang_of_ang(perp_is_rho(), perp, nia); }
	int ipar_of_par(fp par, jtutil::linear_map<fp>::noninteger_action
				nia = jtutil::linear_map<fp>::nia_error)
		{ return my_patch().iang_of_ang(par_is_rho(), par, nia); }

	// ... (perp,par) --> (rho,sigma)
	int irho_of_iperp_ipar(int iperp, int ipar) const
		{ return perp_is_rho() ? iperp : ipar; }
	int isigma_of_iperp_ipar(int iperp, int ipar) const
		{ return perp_is_rho() ? ipar : iperp; }
	fp rho_of_perp_par(fp perp, fp par) const
		{ return perp_is_rho() ? perp : par; }
	fp sigma_of_perp_par(fp perp, fp par) const
		{ return perp_is_rho() ? par : perp; }
	// ... (rho,sigma) --> (perp,par)
	int iperp_of_irho_isigma(int irho, int isigma) const
		{ return perp_is_rho() ? irho : isigma; }
	int ipar_of_irho_isigma(int irho, int isigma) const
		{ return par_is_rho() ? irho : isigma; }
	fp perp_of_rho_sigma(fp rho, fp sigma) const
		{ return perp_is_rho() ? rho : sigma; }
	fp par_of_rho_sigma(fp rho, fp sigma) const
		{ return par_is_rho() ? rho : sigma; }

	// outer perp of nominal grid on this edge
	// ... this is outermost *grid point*
	fp grid_outer_iperp() const
		{ return my_patch().minmax_iang(is_min(), is_rho()); }
	// ... this is actual outer edge of grid
	//     (might be halfway between two grid points)
	fp grid_outer_perp() const
		{ return my_patch().minmax_ang(is_min(), is_rho()); }
	// ... this is grid_outer_perp() converted back to the iperp
	//     coordinate, but still returned as floating-point;
	//     it will be either integer or half-integer
	fp fp_grid_outer_iperp() const
		{ return fp_iperp_of_perp(grid_outer_perp()); }



	//
	// ***** min/max/outer coordinates of edge *****
	//

	// min/max/size ipar of the edge
	// (these are exteme limits for any iperp, a given ghost zone
	//  or interpolation region may have tighter and/or iperp-dependent
	// limits)
	// ... not including corners
	int min_ipar_without_corners() const
		{ return my_patch().min_iang(par_is_rho()); }
	int max_ipar_without_corners() const
		{ return my_patch().max_iang(par_is_rho()); }
	// ... including corners
	int min_ipar_with_corners() const
		{ return my_patch().ghosted_min_iang(par_is_rho()); }
	int max_ipar_with_corners() const
		{ return my_patch().ghosted_max_iang(par_is_rho()); }
	// ... of the corners themselves
	int min_ipar_corner__min_ipar() const
		{ return min_ipar_with_corners(); }
	int min_ipar_corner__max_ipar() const
		{ return min_ipar_without_corners() - 1; }
	int max_ipar_corner__min_ipar() const
		{ return max_ipar_without_corners() + 1; }
	int max_ipar_corner__max_ipar() const
		{ return max_ipar_with_corners(); }

	// membership predicates for ipar corners, non-corners
	bool ipar_is_in_min_ipar_corner(int ipar) const
		{
		return    (ipar >= min_ipar_corner__min_ipar())
		       && (ipar <= min_ipar_corner__max_ipar());
		}
	bool ipar_is_in_max_ipar_corner(int ipar) const
		{
		return    (ipar >= max_ipar_corner__min_ipar())
		       && (ipar <= max_ipar_corner__max_ipar());
		}
	bool ipar_is_in_corner(int ipar) const
		{
		return    ipar_is_in_min_ipar_corner(ipar)
		       || ipar_is_in_max_ipar_corner(ipar);
		}
	bool ipar_is_in_noncorner(int ipar) const
		{
		return    (ipar >= min_ipar_without_corners())
		       && (ipar <= max_ipar_without_corners());
		}

	// convenience function selecting amongst the above
	// membership predicates
	bool ipar_is_in_selected_part(bool want_corners,
				      bool want_noncorner,
				      int ipar)
		const
		{
		return    (want_corners   && ipar_is_in_corner   (ipar))
		       || (want_noncorner && ipar_is_in_noncorner(ipar));
		}

	// outer (farthest from patch center) iperp of nominal grid
	int nominal_grid_outer_iperp() const
		{
		return my_patch()
		       .minmax_iang(is_min(), is_rho());
		}


	//
	// ***** constructor, destructor *****
	//

	patch_edge(patch& my_patch_in,
		   bool is_min_in, bool is_rho_in)
		: my_patch_(my_patch_in),
		  is_min_(is_min_in), is_rho_(is_rho_in)
		{ }
	// compiler-synthesized (no-op) destructor is fine

private:
	// we forbid copying and passing by value
	// by declaring the copy constructor and assignment operator
	// private, but never defining them
	patch_edge(const patch_edge& rhs);
	patch_edge& operator=(const patch_edge& rhs);

private:
	patch& my_patch_;
	const bool is_min_, is_rho_;
	};

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif
