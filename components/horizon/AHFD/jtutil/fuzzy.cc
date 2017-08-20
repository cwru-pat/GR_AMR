// fuzzy.cc -- template for fuzzy comparisons et al on floating point values
// $Header$
//
// jtutil::fuzzy::tolerance - comparison tolerance
// jtutil::fuzzy::EQ
// jtutil::fuzzy::is_integer
// jtutil::fuzzy::floor
// jtutil::fuzzy::ceiling
//
// ***** template instantiations and specializations *****
//

#include <stdlib.h>
#include <cmath>
#include "util.hh"

//******************************************************************************

namespace jtutil
	  {
template <typename fp_t>
bool fuzzy<fp_t>::EQ(fp_t x, fp_t y)
{
fp_t max_abs = jtutil::max(jtutil::abs(x), jtutil::abs(y));
fp_t epsilon = jtutil::max(tolerance_, tolerance_*max_abs);

return jtutil::abs(x-y) <= epsilon;
}
	  }	// namespace jtutil::

//******************************************************************************

namespace jtutil
	  {
template <typename fp_t>
bool fuzzy<fp_t>::is_integer(fp_t x)
{
int i = round<fp_t>::to_integer(x);
return EQ(x, fp_t(i));
}
	  }	// namespace jtutil::

//******************************************************************************

namespace jtutil
	  {
template <typename fp_t>
int fuzzy<fp_t>::floor(fp_t x)
{
return fuzzy<fp_t>::is_integer(x)
       ? round<fp_t>::to_integer(x)
       : round<fp_t>::floor(x);
}
	  }	// namespace jtutil::

//******************************************************************************

namespace jtutil
	  {
template <typename fp_t>
int fuzzy<fp_t>::ceiling(fp_t x)
{
return fuzzy<fp_t>::is_integer(x)
       ? round<fp_t>::to_integer(x)
       : round<fp_t>::ceiling(x);
}
	  }	// namespace jtutil::

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// ***** template instantiations and specializations *****
//

//
// Thanks to Thomas Mang <a9804814@unet.univie.ac.at> for helping
// me figure out the correct syntax here!
//

namespace jtutil
	  {
// initializations of fuzzy::tolerance for each instantation we're going to make
template <>
  float fuzzy<float>::tolerance_ = 1.0e-5;	// about 100 * FLT_EPSILON

template <>
  double fuzzy<double>::tolerance_ = 1.0e-12;	// about 1e4 * DBL_EPSILON

// template instantiations
template class fuzzy<float>;
template class fuzzy<double>;
	  }	// namespace jtutil::
