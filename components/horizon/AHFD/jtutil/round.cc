// round.hh -- template for rounding floating point values
// $Header$
//
// *** Implementation Notes ***
// jtutil::round::to_integer
// jtutil::round::floor
// jtutil::round::ceiling
//
// ***** template instantiations *****
//

#include <stdlib.h>

#include <cmath>
#include "util.hh"

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// *** Implementation Notes ***
//
// We assume throughout this code that C++'s "built-in" conversion
// from <fp_t> to integer takes the floor, at least for zero or positive
// values.
//

//******************************************************************************

// round to nearest integer, up for exact tie
namespace jtutil
	  {
template <typename fp_t>
int round<fp_t>::to_integer(fp_t x)
{
return (x >= 0.0)
       ? int(x + 0.5)		// eg 3.6 --> int(4.1) = 4
       : - int( (-x) + 0.5 );	// eg -3.6 --> - int(4.1) = -4
}
	  }	// namespace jtutil::

//******************************************************************************

namespace jtutil
	  {
template <typename fp_t>
int round<fp_t>::floor(fp_t x)
{
return (x >= 0.0)
       ? int(x)
       : - ceiling(-x);
}
	  }	// namespace jtutil::

//******************************************************************************

namespace jtutil
	  {
template <typename fp_t>
int round<fp_t>::ceiling(fp_t x)
{
return (x >= 0.0)
       ? int(x) + (x != fp_t(int(x)))
       : - floor(-x);
}
	  }	// namespace jtutil::

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// ***** template instantiations *****
//

template class jtutil::round<float>;
template class jtutil::round<double>;
