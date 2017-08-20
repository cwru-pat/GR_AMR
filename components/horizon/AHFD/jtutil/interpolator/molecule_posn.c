 /*@@
   @file      molecule_posn.c
   @date      23 Oct 2001
   @author    Jonathan Thornburg <jthorn@aei.mpg.de>
   @desc      Worker function to compute molecule positions for interpolator.
   @enddesc
   @version   $Id$
 @@*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "../../AHFD_macros.h"

#include "InterpLocalUniform.h"

/* the rcs ID and its dummy function to use it */
/* static const char *rcsid = "$Header$"; */
/* #ifndef AEILOCALINTERP_STANDALONE_TEST */
/*   CCTK_FILEVERSION(AEIThorns_AEILocalInterp_src_molecule_posn_c) */
/* #endif */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*@@
 @routine	AEILocalInterp_molecule_posn
 @date		23 Oct 2001
 @author	Jonathan Thornburg <jthorn@aei.mpg.de>
 @desc
	Given a uniformly-spaced grid, this function computes molecule
	positions for an interpolation or similar operation:
	(1) It computes the closest grid point to the input coordinate
	(2) It does the slightly tricky odd/even computation to decide where
	    to center the molecule
	(3) It checks for the molecule falling off the edge of the grid,
	    and off-centers it as appropriate.

	For (1), we assume that grid points have floating-point
	coordinates (we refer to these as "x") which are an arbitrary
	linear function of the integer grid coordinates (we refer to
	these as "i"), x = grid_origin + i*grid_delta.

	For (2), suppose we have a data array indexed by [i], from which we
	wish to select an N-point molecule [i_lo, i_hi] centered as close as
	possible to the point  x .  The problem is just how to choose the
	centering.  The following diagram illustrates the issues:

	N  i_lo  i_hi  [i-3]  [i-2]  [i-1]   [i]   [i+1]  [i+2]  [i+3]  [i+4]
	-  ----  ----  -----  -----  -----  -----  -----  -----  -----  -----
	2  i     i+1                          *xxxxxx*
	3  i-1   i+1                   *   xxx*xxx   *
	4  i-1   i+2                   *      *xxxxxx*      *
	5  i-2   i+2            *      *   xxx*xxx   *      *
	6  i-2   i+3            *      *      *xxxxxx*      *      *
	7  i-3   i+3     *      *      *   xxx*xxx   *      *      *
	8  i-3   i+4     *      *      *      *xxxxxx*      *      *      *

	The diagram shows the range of  x  values relative to the grid,
	for which we'll choose each centering.  Thus if N is even, the
	centering is based on which grid zone contains x (take the floor
	of the floating-point i coordinate), while if N is odd, the centering
	is based on which grid point is closest to x (round the floating-
	-point i coordinate to the nearest integer).

	It's also convenient to introduce the integer molecule coordinate m;
	this is the integer grid coordinate i relative to that of the molecule
	center, i.e. the above diagram has columns labelled with [i+m].

	For (3) above, note that we describe the range of the grid
	by the *closed* interval [grid_i_min, grid_i_max].
 @enddesc

 @hdate    27.Jan.2003
 @hauthor  Jonathan Thornburg <jthorn@aei.mpg.de>
 @hdesc    Complete rewrite: now supports
             @var boundary_off_centering_tolerance_min,
             @var boundary_off_centering_tolerance_max,
             @var boundary_extrapolation_tolerance_min,
             @var boundary_extrapolation_tolerance_max,
           also change to returning status code,
           also drop returning @var min_m and @var max_m
           since they were never used.
 @endhdesc 

 @var	grid_origin
 @vdesc	The floating-point coordinate x of the grid point i=0.
 @vtype	fp grid_origin
 @endvar

 @var	grid_delta
 @vdesc	The grid spacing (in the floating-point coordinate x).
 @vtype	fp grid_delta
 @endvar

 @var	grid_i_min
 @vdesc	The minimum integer coordinate i of the grid.
 @vtype	int grid_i_min
 @endvar

 @var	grid_i_max
 @vdesc	The maximum integer coordinate i of the grid.
 @vtype	int grid_i_max
 @endvar

 @var	molecule_size
 @vdesc	The size (number of points) of the molecule.
 @vtype	int molecule_size
 @endvar

 @var	boundary_off_centering_tolerance_min,
	boundary_off_centering_tolerance_max
 @vdesc	Specifies how many grid spacings the interpolation point
	may be beyond the default-centering region before being
	declared "out of range" on the {minimum,maximum} boundaries
	of the grid respectively.
 @vtype fp boundary_off_centering_tolerance_min,
	   boundary_off_centering_tolerance_max
 @endvar

 @var	boundary_extrapolation_tolerance_min,
	boundary_extrapolation_tolerance_max
 @vdesc	Specifies how many grid spacings the interpolation point
	may be beyond the grid boundary before being declared
	"out of range" on the {minimum,maximum} boundaries
	of the grid respectively.
 @vtype fp boundary_extrapolation_tolerance_min,
	   boundary_extrapolation_tolerance_max
 @endvar

 @var	x
 @vdesc	The floating-point coordinate of the interpolation point.
 @vtype	fp x
 @endvar

 @var	debug
 @vdesc	A debugging flag (0 = no debug output, > 0 = print debug output)
 @vtype	int x
 @endvar

 @var	i_center
 @vdesc	A pointer to an value where this function should
	store the integer coordinate of the molecule center,
	or NULL to skip storing this
 @vtype	int *i_center
 @vio	pointer to out
 @endvar

 @var	x_rel
 @vdesc	A pointer to where this function should store the
	interpolation point's x coordinate relative to the
	molecule center, measured in units of the grid spacing;
	or NULL to skip storing this.
 @vtype	fp *x_rel
 @vio	pointer to out
 @endvar

 @returntype	int
 @returndesc
	This function returns 0 if the interpolation point is "in range",
	or one of the (negative) error codes defined in "InterpLocalUniform.h"
	if the interpolation point is "out of range":
	MOLECULE_POSN_ERROR_X_LT_MIN
		if  x  is out-of-range on the min end of the grid
		(i.e. x < the minimum allowable x)
	MOLECULE_POSN_ERROR_X_GT_MAX
		if  x  is out-of-range on the max end of the grid
		(i.e. x > the maximum allowable x)
	MOLECULE_POSN_ERROR_GRID_TINY
		if the grid is smaller than the molecule
	MOLECULE_POSN_ERROR_NAN
		if we encounter a NaN or other non-finite floating-point value
	MOLECULE_POSN_ERROR_DX_ZERO
		if  grid_delta == 0.0
 @endreturndesc
 @@*/
int AEILocalInterp_molecule_posn(fp grid_origin, fp grid_delta,
				 int grid_i_min, int grid_i_max,
				 int molecule_size,
				 fp boundary_off_centering_tolerance_min,
				 fp boundary_off_centering_tolerance_max,
				 fp boundary_extrapolation_tolerance_min,
				 fp boundary_extrapolation_tolerance_max,
				 fp x,
				 int debug,
				 int* i_center, fp* x_rel)
{
/*
 * ***** IMPORTANT *****
 *
 * This code is ++tricky.  (More accurately, the basic algorithm is
 * fairly simple, but there are lots of corner cases.)  There's a
 * fairly rigorous test driver for it in the file  test_molecule_posn.c
 * in this directory.  This is a standalone test driver, which can be
 * compiled via
 *	gmake -f Makefile.standalone
 * in this directory.  If you make any changes to this code, please
 * verify that all the tests still pass!
 */

if (debug >= 8)
   then {
	printf("AEILocalInterp_molecule_posn():\n");
	printf("   grid_origin=%g grid_delta=%g\n",
	       (double) grid_origin, (double) grid_delta);
	printf("   grid_i_[min,max]=[%d,%d] molecule_size=%d\n",
	       grid_i_min, grid_i_max, molecule_size);
	printf("   boundary_off_centering_tolerance_[min,max]=[%g,%g]\n",
	       (double) boundary_off_centering_tolerance_min,
	       (double) boundary_off_centering_tolerance_max);
	printf("   boundary_extrapolation_tolerance_[min,max]=[%g,%g]\n",
	       (double) boundary_extrapolation_tolerance_min,
	       (double) boundary_extrapolation_tolerance_max);
	printf("   x=%g\n", (double) x);
	}

/*
 * basic sanity checks
 */

/* is the molecule larger than the grid? */
if (molecule_size > HOW_MANY_IN_RANGE(grid_i_min,grid_i_max))
   then return MOLECULE_POSN_ERROR_GRID_TINY;		/*** ERROR RETURN ***/

/* has someone passed us NaNs for coordinates etc? */
if (!isfinite(grid_origin) || !isfinite(grid_delta) || !isfinite(x))
   then return MOLECULE_POSN_ERROR_NAN;
if (    !isfinite(boundary_off_centering_tolerance_min)
     || !isfinite(boundary_off_centering_tolerance_max)
     || !isfinite(boundary_extrapolation_tolerance_min)
     || !isfinite(boundary_extrapolation_tolerance_max)    )
   then return MOLECULE_POSN_ERROR_NAN;

/* is the grid spacing zero (we'll need to divide by it)? */
if (grid_delta == 0.0)
   then return MOLECULE_POSN_ERROR_DX_ZERO;


/* molecule radia (inherently positive) in +/- directions */
const int mr_plus  = (molecule_size >> 1);
const int mr_minus = molecule_size - mr_plus - 1;

/* range of x_rel for which this molecule is centered, */
/* cf. diagram in header comment */
const fp centered_min_x_rel = IS_EVEN(molecule_size) ? 0.0 : -0.5;
const fp centered_max_x_rel = IS_EVEN(molecule_size) ? 1.0 : +0.5;

/* range of i where a centered molecule would fit within the grid, */
/* as floating-point numbers */
const fp fp_centered_min_possible_i = grid_i_min + mr_minus
						 + centered_min_x_rel;
const fp fp_centered_max_possible_i = grid_i_max - mr_plus
						 + centered_max_x_rel;

/* integer coordinate i of interpolation point, as a floating-point number */
const fp fp_i = (x - grid_origin) / grid_delta;
assert(isfinite(fp_i));

if (debug > 8)
   then {
	printf("   mr_{plus,minus}={%d,%d}\n", mr_plus, mr_minus);
	printf("   centered_[min,max]_x_rel=[%g,%g]\n",
	       (double) centered_min_x_rel, (double) centered_max_x_rel);
	printf("   fp_centered_[min,max]_possible_i=[%g,%g]\n",
	       (double) fp_centered_min_possible_i,
	       (double) fp_centered_max_possible_i);
	printf("   fp_i=%g\n", (double) fp_i);
	}

/* is interpolation point x beyond the extrapolation tolerance? */
if (fp_i < (fp)grid_i_min - boundary_extrapolation_tolerance_min)
   then return MOLECULE_POSN_ERROR_X_LT_MIN;		/*** ERROR RETURN ***/
if (fp_i > (fp)grid_i_max + boundary_extrapolation_tolerance_max)
   then return MOLECULE_POSN_ERROR_X_GT_MAX;		/*** ERROR RETURN ***/

/* is interpolation point x beyond the off-centering tolerance? */
if (fp_i < fp_centered_min_possible_i - boundary_off_centering_tolerance_min)
   then return MOLECULE_POSN_ERROR_X_LT_MIN;		/*** ERROR RETURN ***/
if (fp_i > fp_centered_max_possible_i + boundary_off_centering_tolerance_max)
   then return MOLECULE_POSN_ERROR_X_GT_MAX;		/*** ERROR RETURN ***/

/* now choose the actual molecule position/centering: */
/* first set up a centered molecule */
  {
/* compute as a floating point number */
/* (because that's what the C math library provides), */
/* but this value is guaranteed to be integral */
fp fp_i_center = IS_EVEN(molecule_size)
		 ? floor(fp_i)
		 : ROUND_TO_INTEGER__RESULT_IS_FP(fp_i);
assert(isfinite(fp_i_center));
int int_i_center = (int) fp_i_center;	/* convert to integer */
if (debug > 8)
   then printf("   initial: fp_i_center=%g int_i_center=%d\n",
	       (double) fp_i_center, int_i_center);

/* clamp the molecule at the grid boundaries */
if (int_i_center < grid_i_min + mr_minus)
   then {
	int_i_center = grid_i_min + mr_minus;
	fp_i_center = (fp) int_i_center;
	}
if (int_i_center > grid_i_max - mr_plus)
   then {
	int_i_center = grid_i_max - mr_plus;
	fp_i_center = (fp) int_i_center;
	}

if (debug > 8)
   then printf(
	"   final result after clamping: fp_i_center=%g int_i_center=%d\n",
	       (double) fp_i_center, int_i_center);

/* store the results */
if (i_center != NULL)
   then *i_center = int_i_center;
if (x_rel != NULL)
   then *x_rel = fp_i - fp_i_center;

return 0;						/*** NORMAL RETURN ***/
  }
}
