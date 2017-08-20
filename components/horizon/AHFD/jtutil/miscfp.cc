// miscfp.cc -- misc floating-point functions

//
// jtutil::signum - sign of a floating point number
// jtutil::hypot3 - 3D Euclidean distance
// jtutil::arctan_xy - 4-quadrant arc tangent
// jtutil::modulo_reduce - reduce an angle modulo 2*pi radians (360 degrees)
// jtutil::zero_C_array - set a C-style array to all zeros
//
// ***** template instantiations *****
//

#include <math.h>
#include <stdlib.h>

// we want to instantiate templates with CCTK_* types
#include "../cctk.h"

#include "util.hh"
using namespace SAMRAI;
//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function computes the floating point "signum" function (as in APL),
// signum(x) = -1.0, 0.0, or +1.0, according to the sign of x.
//
namespace jtutil
	  {
double signum(double x)
{
if (x == 0.0)
   then return 0.0;
   else return (x > 0.0) ? 1.0 : -1.0;
}
	  }	// namespace jtutil::

//******************************************************************************

//
// This function computes the 3-dimensional Euclidean distance,
// analogously to the standard math library function  hypot(2) .
//
// Arguments:
// (x,y,z) = (in) The rectangular coordinates of a point in $\Re^3$.
//
// Result:
// The function returns the Euclidean distance of (x,y,z) from the origin.
//
// Bugs:
// Unlike  hypot(2), this function takes no special care to avoid
// unwarranted IEEE exceptions if any of |x|, |y|, or |z| is close to
// the overflow and/or underflow threshold.
//
namespace jtutil
	  {
double hypot3(double x, double y, double z)
{
return sqrt(x*x + y*y + z*z);
}
	  }	// namespace jtutil::

//******************************************************************************

//
// This function computes a four-quadrant arctangent, using what I think
// is the "right" ordering of the arguments, and returning 0.0 if both
// arguments are 0.0.
//
// Arguments:
// (x,y) = (in) The rectangular coordinates of a point in $\Re^2$.
//
// Result:
// The function returns a value $\theta$ such that $x + iy = R e^{i\theta}$ for
// some real $R$, i.e. it returns the angle between the positive $x$ axis and
// the line joining the origin and the point $(x,y)$.
//
namespace jtutil
	  {
double arctan_xy(double x, double y)
{
// note reversed argument order (y,x) in std::atan2() function
return ((x == 0.0) && (y == 0.0)) ? 0.0 : atan2(y,x);
}
	  }	// namespace jtutil::

//******************************************************************************

//
// This function reduces  x  modulo  xmod  to be (fuzzily) in the range
//  [xmin, xmax] , or does an  error_exit()  if no such value exists.
//
namespace jtutil
	  {
double modulo_reduce(double x, double xmod, double xmin, double xmax)
{
double xx = x;

	while (fuzzy<double>::LT(xx, xmin))
	{
	xx += xmod;
	}

	while (fuzzy<double>::GT(xx, xmax))
	{
	xx -= xmod;
	}

if (! (fuzzy<double>::GE(xx, xmin) && fuzzy<double>::LE(xx, xmax)) )
   error_exit(ERROR_EXIT,
"***** modulo_reduce(): no modulo value is fuzzily within specified range!\n"
"                       x = %g   xmod = %g\n"
"                       [xmin,xmax] = [%g,%g]\n"
"                       ==> xx = %g\n"
,
			   x, xmod,
			   xmin, xmax,
			   xx);					/*NOTREACHED*/

return xx;
}
	  }	// namespace jtutil::

//******************************************************************************

//
// This function sets a C-style array to all zeros.
//
namespace jtutil
	  {
template <typename fp_t>
  void zero_C_array(int N, fp_t array[])
{
	for (int i = 0 ; i < N ; ++i)
	{
	array[i] = 0;
	}
}
	  }	// namespace jtutil::

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// ***** template instantiations *****
//

namespace jtutil
	  {
template void zero_C_array<CCTK_REAL>(int, CCTK_REAL[]);
	  }	// namespace jtutil::
