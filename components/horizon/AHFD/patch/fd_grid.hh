#ifndef AHFD_PATCH_FD_GRID_H
#define AHFD_PATCH_FD_GRID_H

// fd_grid.hh -- grid with finite differencing operations
// $Header$
//
#include <cstdio>
#include <assert.h>
#include <math.h>
// #include "../jtutil/util.hh"
// #include "../jtutil/array.hh"
// #include "../jtutil/linear_map.hh"
//#include "coords.hh"
#include "grid.hh"
// *** Implementation Notes -- Overview ***
// *** Implementation Notes -- Techniques using C++ Templates ***
// *** Implementation Notes -- Techniques using the C/C++ Preprocessor ***
// *** Implementation Notes -- Run-Time Choice of Molecules ***
// *** finite difference molecules ***
// ***** fd_grid - grid with finite differencing operations *****
//

//
// prerequisites:
//    <stdio.h>
//    <assert.h>
//    <math.h>
//    "stdc.h"
//    "config.hh"
//    "../jtutil/util.hh"
//    "../jtutil/array.hh"
//    "../jtutil/linear_map.hh"
//    "coords.hh"
//    "grid.hh"
//

// everything in this file is inside this namespace
namespace AHFD
	  {

//******************************************************************************

//
// *** Implementation Notes -- Overview ***
//

//
// The key design problem for our finite differencing is how to
// implement an entire family of 5(9) finite difference operations in
// 2D(3D)
//
//	partial_rho		partial_sigma
//	partial_{rho,rho}	partial_{rho,sigma}
//				partial_{sigma,sigma}
//
//	partial_x		partial_y		partial_z
//	partial_xx		partial_xy		partial_xz
//				partial_yy		partial_yz
//							partial_zz
//
// without having to write out the finite differencing molecules multiple
// times, and while still preserving maximum inline-function efficiency.
// In particular, mixed 2nd-order derivative operations like partial_xy
// should be automatically composed from the two individual 1st derivative
// operations (partial_x and partial_y).
//

//
// Our basic approach is to define each finite difference molecule in
// a generic 1-dimensional form using an abstract "data(m)" interface.
// Here we use the terminology that a finite difference molecule is
// defined as
//	out[k] = sum(m) c[m] * in[k+m]
// where c[] is the vector/matrix of molecule coefficients, and m is
// the (integer) relative grid coordinate within a molecule.
//
// That is, for example, we define the usual 2nd order centered 1st
// derivative operator as
//	diff = 0.5*inv_delta_x*(data(+1) - data(-1))
// leaving unspecified just what the data source is.  We then use this
// with an appropriate data source (indexing along that gridfn array axis)
// for each directional derivative operation, and we compose two of
// these, using the first along x as the data source for the second
// along y, for the mixed 2nd-order derivative operation.
//

//******************************************************************************

//
// *** Implementation Notes -- Techniques using C++ Templates ***
// 

//
// There are two plausible ways to use C++ templates
//	[C++ templates are described in detail in chapter 13 of
//	Stroustrup's "The C++ Programming Language" (3rd Edition),
//	hereinafter "C++PL", and chapter 15 of Stroustrup's
//	"The Design and Evolution of C++", hereinafter "D&EC++".]
// to write the sort of generic-at-compile-time code we want:
// - Template specializations for each axis, as discussed in D&EC++
//   section 15.10.3.
// - Overloaded functions for each axis, with an argument type
//   (possibly that of an extra unused argument) selecting the
//   appropriate axis and hence the appropriate function.  This
//   technique is discussed in D&EC++ section 15.6.3.1.
//
// Quoting from D&EC++ (section 15.6.3.1),
//
//	The fundamental observation is that every property
//	of a type or an algorithm can be represented by a
//	type (possibly defined specificaly to do exactly
//	that).  That done, such a type can be used to guide
//	the overload resolution to select a function that
//	depends on the desired property.  [...]
//
//	Please note that thanks to inlining this resolution
//	is done at compile-time, so the appropriate [...]
//	function will be called directly without any run-time
//	overhead.
//
// Quoting from C++PL3 (section 13.4),
//
//	Passing [...] operations as a template parameter has two
//	significant benefits compared to alternatives such as
//	passing pointers to functions.  Several operations can
//	be passed as a single argument with no run-time cost.
//	In addition, the [...] operators [passed this way] are
//	trivial to inline, whereas inlininkg a call through a
//	pointer to function requires exceptional attention from
//	a compiler.
//

//
// In my opinion the template-specialization design is cleaner, and it
// clearly has no run-time cost (whereas the overloaded-function design
// may have a run-time cost for constructing and passing unused objects),
// so we use it here.
//
// There are, however, two (non-fatal) problema with this approach:
// - Unfortunately, it appears C++ (or at least gcc 2.95.1) forbids
//   template specialization within a class, so some of the functions
//   which whould logically be class members, must instead be defined
//   outside any class.  We use the namespace  fd_stuff::  to hide
//   these from the outside world.
// - C++PL3, section C.13.3, states that
//	Only class templates can be template arguments.
//   so we have to use dummy classes around some of our template
//   functions.  To avoid extra constructor/destructor overhead, we
//   make these template functions static.
//

//******************************************************************************

//
// *** Implementation Notes -- Techniques using the C/C++ Preprocessor ***
//

//
// The fundamental problem with the template approaches is portability:
// Although the C++ standard describes powerful template facilities, not
// all C++ compilers yet fully support these.  As an alternative, we can
// use the C/C++ preprocessor.  This is ugly and dangerous (global names!),
// but is probably simpler than any of the template approaches.  It can
// provide the same finite differencing functionality and efficiency as
// the template-based approaches.
//
// Because of its greater portability, we use the preprocessor-based
// approach here.
//

//******************************************************************************


//
// *** Implementation Notes -- Run-Time Choice of Molecules ***
//
// *If* we want to allow the finite differencing scheme to be changed
// at run-time (e.g. from a parameter file), there are three plausible
// ways to do this:
// - Using  switch(molecule_type) , as is standard in C.  This is
//   simple, and for this particular application quite well-structured
//   and maintainable (there are only a few different molecule types,
//   all centralized in this file).
// - Using virtual functions, with  molecule  a virtual base class
//   and individual molecules derived from it.  This is elegant, but
//   may have some performance problems (below).  It also requires some
//   sort of switch-based "object factory" to interface with with the
//   molecule-choice parameters.
// - Write all the finite differencing code multiple times, once for
//   each finite differencing scheme.
//
// The typical use of these functions will be from within a loop over
// a whole grid.  In both cases we can expect excellent accuracy from
// modern hardware branch prediction (and thus minimal performance loss
// from the branching).  It's reasonable to expect a compiler to fully
// inline the switch-based code, exposing all the gridfn array subscriptings
// to strength reduction etc, but this is much trickier for the
// virtual-function--based code.  For this reason, the switch-based
// design seems superior to the virtual-function--based one.
//
// However, at present we don't implement any run-time selection: we
// "just" fix the finite differencing scheme at compile time via the
// preprocessor.
//

//******************************************************************************

//
// *** finite difference molecules ***
//

//**************************************

//
// define the actual molecules
//
// In the following macros, we first define all the distinct floating-
// -point numbers appearing in a molecules as "K" constants (all > 0),
// then define the actual derivative and its molecule coefficients
// using +/- the "K" constants, with multiplies by 1.0 elided and 0
// terms skipped in computing the derivative.  This (hopefully) gives
// maximum efficiency by avoiding the generated code loading the same
// constants multiple times.
// 

//
// The molecule macros all take the following arguments:
// inv_delta_x_ = inverse of grid spacing in the finite differencing
//		  direction
// data_= a data-fetching function or macro: data_(ghosted_gfn, irho, isigma)
//	  is the data to be finite differenced
// irho_plus_m_ = a function or macro: irho_plus_m_(irho,m) returns the
//		  rho coordinate to be passed to data_() for the [m]
//		  molecule coefficient
// isigma_plus_m_ = same thing, for the sigma coordinate
//
// n.b. We grab the variables ghosted_gfn, irho, and isigma from the calling
//      environment, and we define assorted local variables as needed!
//

//**************************************

//
// 2nd order
//

#define FD_GRID__ORDER2__MOL_RADIUS	1
#define FD_GRID__ORDER2__MOL_DIAMETER	3

#define FD_GRID__ORDER2__DX__KPM1	0.5
#define FD_GRID__ORDER2__DX(inv_delta_x_, data_,			\
			    irho_plus_m_, isigma_plus_m_)		\
	const fp data_p1 = data_(ghosted_gfn,				\
				   irho_plus_m_(irho  ,+1),		\
				 isigma_plus_m_(isigma,+1));		\
	const fp data_m1 = data_(ghosted_gfn,				\
				   irho_plus_m_(irho  ,-1),		\
				 isigma_plus_m_(isigma,-1));		\
	const fp sum = FD_GRID__ORDER2__DX__KPM1 * (data_p1 - data_m1);	\
	return inv_delta_x_ * sum;				/* end macro */
#define FD_GRID__ORDER2__DX__COEFF_M1	(-FD_GRID__ORDER2__DX__KPM1)
#define FD_GRID__ORDER2__DX__COEFF_0	0.0
#define FD_GRID__ORDER2__DX__COEFF_P1	(+FD_GRID__ORDER2__DX__KPM1)

#define FD_GRID__ORDER2__DXX__K0		2.0
#define FD_GRID__ORDER2__DXX(inv_delta_x_, data_,			\
			     irho_plus_m_, isigma_plus_m_)		\
	const fp data_p1 = data_(ghosted_gfn,				\
				   irho_plus_m_(irho  ,+1),		\
				 isigma_plus_m_(isigma,+1));		\
	const fp data_0  = data_(ghosted_gfn,				\
				   irho_plus_m_(irho  , 0),		\
				 isigma_plus_m_(isigma, 0));		\
	const fp data_m1 = data_(ghosted_gfn,				\
				   irho_plus_m_(irho  ,-1),		\
				 isigma_plus_m_(isigma,-1));		\
	const fp sum							\
		= data_m1 - FD_GRID__ORDER2__DXX__K0*data_0 + data_p1;	\
	return jtutil::pow2(inv_delta_x_) * sum;		/* end macro */
#define FD_GRID__ORDER2__DXX__COEFF_M1	1.0
#define FD_GRID__ORDER2__DXX__COEFF_0	(-FD_GRID__ORDER2__DXX__K0)
#define FD_GRID__ORDER2__DXX__COEFF_P1	1.0

//**************************************

//
// 4th order
//

#define FD_GRID__ORDER4__MOL_RADIUS	2
#define FD_GRID__ORDER4__MOL_DIAMETER	5

#define FD_GRID__ORDER4__DX__KPM2	(1.0/12.0)
#define FD_GRID__ORDER4__DX__KPM1	(8.0/12.0)
#define FD_GRID__ORDER4__DX(inv_delta_x_, data_,			\
			    irho_plus_m_, isigma_plus_m_)		\
	const fp data_p2 = data_(ghosted_gfn,				\
				   irho_plus_m_(irho  ,+2),		\
				 isigma_plus_m_(isigma,+2));		\
	const fp data_p1 = data_(ghosted_gfn,				\
				   irho_plus_m_(irho  ,+1),		\
				 isigma_plus_m_(isigma,+1));		\
	const fp data_m1 = data_(ghosted_gfn,				\
				   irho_plus_m_(irho  ,-1),		\
				 isigma_plus_m_(isigma,-1));		\
	const fp data_m2 = data_(ghosted_gfn,				\
				   irho_plus_m_(irho  ,-2),		\
				 isigma_plus_m_(isigma,-2));		\
	const fp sum							\
		=   FD_GRID__ORDER4__DX__KPM1 * (data_p1 - data_m1)	\
		  + FD_GRID__ORDER4__DX__KPM2 * (data_m2 - data_p2);	\
	return inv_delta_x_ * sum;				/* end macro */
#define FD_GRID__ORDER4__DX__COEFF_M2	(+FD_GRID__ORDER4__DX__KPM2)
#define FD_GRID__ORDER4__DX__COEFF_M1	(-FD_GRID__ORDER4__DX__KPM1)
#define FD_GRID__ORDER4__DX__COEFF_0	0.0
#define FD_GRID__ORDER4__DX__COEFF_P1	(+FD_GRID__ORDER4__DX__KPM1)
#define FD_GRID__ORDER4__DX__COEFF_P2	(-FD_GRID__ORDER4__DX__KPM2)

//**************************************

#define FD_GRID__ORDER4__DXX__KPM2	( 1.0/12.0)
#define FD_GRID__ORDER4__DXX__KPM1	(16.0/12.0)
#define FD_GRID__ORDER4__DXX__K0	(30.0/12.0)
#define FD_GRID__ORDER4__DXX(inv_delta_x_, data_,			\
			     irho_plus_m_, isigma_plus_m_)		\
	const fp data_p2 = data_(ghosted_gfn,				\
				   irho_plus_m_(irho  ,+2),		\
				 isigma_plus_m_(isigma,+2));		\
	const fp data_p1 = data_(ghosted_gfn,				\
				   irho_plus_m_(irho  ,+1),		\
				 isigma_plus_m_(isigma,+1));		\
	const fp data_0  = data_(ghosted_gfn,				\
				   irho_plus_m_(irho  , 0),		\
				 isigma_plus_m_(isigma, 0));		\
	const fp data_m1 = data_(ghosted_gfn,				\
				   irho_plus_m_(irho  ,-1),		\
				 isigma_plus_m_(isigma,-1));		\
	const fp data_m2 = data_(ghosted_gfn,				\
				   irho_plus_m_(irho  ,-2),		\
				 isigma_plus_m_(isigma,-2));		\
	const fp sum							\
		= - FD_GRID__ORDER4__DXX__K0 * data_0			\
		  + FD_GRID__ORDER4__DXX__KPM1 * (data_m1 + data_p1)	\
		  - FD_GRID__ORDER4__DXX__KPM2 * (data_m2 + data_p2);	\
	return jtutil::pow2(inv_delta_x_) * sum;		/* end macro */
#define FD_GRID__ORDER4__DXX__COEFF_M2	(-FD_GRID__ORDER4__DXX__KPM2)
#define FD_GRID__ORDER4__DXX__COEFF_M1	(+FD_GRID__ORDER4__DXX__KPM1)
#define FD_GRID__ORDER4__DXX__COEFF_0	(-FD_GRID__ORDER4__DXX__K0  )
#define FD_GRID__ORDER4__DXX__COEFF_P1	(+FD_GRID__ORDER4__DXX__KPM1)
#define FD_GRID__ORDER4__DXX__COEFF_P2	(-FD_GRID__ORDER4__DXX__KPM2)

//******************************************************************************

//
// choose finite differencing order via preprocessor symbol FINITE_DIFF_ORDER
//
#ifndef FINITE_DIFF_ORDER
  #error "must define FINITE_DIFF_ORDER!"
#endif

#if   (FINITE_DIFF_ORDER == 2)
  #define FD_GRID__MOL_RADIUS	FD_GRID__ORDER2__MOL_RADIUS
  #define FD_GRID__MOL_DIAMETER	FD_GRID__ORDER2__MOL_DIAMETER
  #define FD_GRID__DX		FD_GRID__ORDER2__DX
  #define FD_GRID__DXX		FD_GRID__ORDER2__DXX
#elif (FINITE_DIFF_ORDER == 4)
  #define FD_GRID__MOL_RADIUS	FD_GRID__ORDER4__MOL_RADIUS
  #define FD_GRID__MOL_DIAMETER	FD_GRID__ORDER4__MOL_DIAMETER
  #define FD_GRID__DX		FD_GRID__ORDER4__DX
  #define FD_GRID__DXX		FD_GRID__ORDER4__DXX
#else
  #error "unknown value " FINITE_DIFF_ORDER " for FINITE_DIFF_ORDER!"
#endif

#define FD_GRID__MOL_AREA	(FD_GRID__MOL_DIAMETER * FD_GRID__MOL_DIAMETER)

//******************************************************************************

//
// ***** fd_grid - grid with finite differencing operations *****
//
// An  fd_grid  is identical to a  grid  except that it also defines
// (rho,sigma)-coordinate finite differencing operations on gridfns.
//

class	fd_grid
	: public grid
	{
	//
	// molecule sizes
	//
public:
	// n.b. this interface implicitly assumes that all molecules
	//      are centered and are the same order and size
	static
	  int finite_diff_order() { return FINITE_DIFF_ORDER; }
	static
	  int molecule_radius()   { return FD_GRID__MOL_RADIUS; }
	static
	  int molecule_diameter() { return FD_GRID__MOL_DIAMETER; }
	static
	  int molecule_min_m() { return -FD_GRID__MOL_RADIUS; }
	static
	  int molecule_max_m() { return  FD_GRID__MOL_RADIUS; }

	//
	// helper functions to compute (irho,isigma) + [m]
	// along each axis
	//
private:
	static
	  int     rho_axis__irho_plus_m(int irho  , int m) { return irho  +m; }
	static
	  int   rho_axis__isigma_plus_m(int isigma, int m) { return isigma  ; }
	static
	  int   sigma_axis__irho_plus_m(int irho  , int m) { return irho    ; }
	static
	  int sigma_axis__isigma_plus_m(int isigma, int m) { return isigma+m; }


	//
	// ***** finite differencing *****
	//
public:

	// 1st derivatives
	fp partial_rho(int ghosted_gfn,  int irho, int isigma)
		const
		{
		FD_GRID__DX(inverse_delta_rho(),
			    ghosted_gridfn,
			    rho_axis__irho_plus_m,
			    rho_axis__isigma_plus_m);
		}
	fp partial_sigma(int ghosted_gfn,  int irho, int isigma)
		const
		{
		FD_GRID__DX(inverse_delta_sigma(),
			    ghosted_gridfn,
			    sigma_axis__irho_plus_m,
			    sigma_axis__isigma_plus_m);
		}

	// "pure" 2nd derivatives
	fp partial_rho_rho(int ghosted_gfn,  int irho, int isigma)
		const
		{
		FD_GRID__DXX(inverse_delta_rho(),
			     ghosted_gridfn,
			     rho_axis__irho_plus_m,
			     rho_axis__isigma_plus_m);
		}
	fp partial_sigma_sigma(int ghosted_gfn,  int irho, int isigma)
		const
		{
		FD_GRID__DXX(inverse_delta_sigma(),
			     ghosted_gridfn,
			     sigma_axis__irho_plus_m,
			     sigma_axis__isigma_plus_m);
		}

	// mixed 2nd partial derivative
	fp partial_rho_sigma(int ghosted_gfn,   int irho, int isigma)
		const
		{
		FD_GRID__DX(inverse_delta_rho(),
			    partial_sigma,
			    rho_axis__irho_plus_m,
			    rho_axis__isigma_plus_m);
		}


	//
	// ***** molecule coefficients *****
	//
public:
	// molecule coefficients
	// n.b. this interface implicitly assumes that all molecules
	//      are position-independent
	fp partial_rho_coeff        (int m) const
		{ return inverse_delta_rho()   * dx_coeff(m); }
	fp partial_sigma_coeff      (int m) const
		{ return inverse_delta_sigma() * dx_coeff(m); }
	fp partial_rho_rho_coeff    (int m) const
		{ return jtutil::pow2(inverse_delta_rho()) * dxx_coeff(m); }
	fp partial_sigma_sigma_coeff(int m) const
		{ return jtutil::pow2(inverse_delta_sigma()) * dxx_coeff(m); }
	fp partial_rho_sigma_coeff(int m_rho, int m_sigma) const
		{
		return partial_rho_coeff(m_rho) * partial_sigma_coeff(m_sigma);
		}

	// worker functions: molecule coefficients for unit grid spacing
private:
	static
	  fp dx_coeff(int m);
	static
	  fp dxx_coeff(int m);


	//
	// ***** constructor, destructor *****
	//
public:

	// constructor: pass through to grid:: constructor
	fd_grid(const grid_array_pars& grid_array_pars_in,
		const grid_pars& grid_pars_in)
		: grid(grid_array_pars_in, grid_pars_in)
		{ }
	// compiler-generated default destructor is ok


private:
        // we forbid copying and passing by value
        // by declaring the copy constructor and assignment operator
        // private, but never defining them
	fd_grid(const fd_grid& rhs);
	fd_grid& operator=(const fd_grid& rhs);
	};

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif
