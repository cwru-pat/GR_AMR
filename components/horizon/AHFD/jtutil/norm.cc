// norm.cc -- compute norms
// $Header$

//
// jtutil::norm::data
// jtutil::norm::mean
// jtutil::norm::two_norm
// jtutil::norm::rms_norm
//
// ***** template instantiations *****
//

#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include "util.hh"

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

//
// This function constructs a  norm  object.
//
namespace jtutil
	  {
template <typename fp_t>
  norm<fp_t>::norm()
	: N_(0L),
	  sum_(0.0), sum2_(0.0),
	  max_abs_value_(0.0), min_abs_value_(0.0),
	  max_value_(0.0), min_value_(0.0)
{ }
	  }

//*****************************************************************************

//
// This function resets a  norm  object to its initial state
//
namespace jtutil
	  {
template <typename fp_t>
  void norm<fp_t>::reset()
{
N_ = 0L;
sum_ = 0.0;
sum2_ = 0.0;
max_abs_value_ = 0.0;
min_abs_value_ = 0.0;
max_value_ = 0.0;
min_value_ = 0.0;
}
	  }

//*****************************************************************************

//
// This function updates the norms with a new data point  x .
//
namespace jtutil
	  {
template <typename fp_t>
  void norm<fp_t>::data(fp_t x)
{
sum_  += x;
sum2_ += x*x;

const fp_t abs_x = jtutil::abs<fp_t>(x);
max_abs_value_ = jtutil::max(max_abs_value_, abs_x);
min_abs_value_ = (N_ == 0) ? abs_x : jtutil::min(min_abs_value_, abs_x);

min_value_ = (N_ == 0) ? x : jtutil::min(min_value_, x);
max_value_ = (N_ == 0) ? x : jtutil::max(max_value_, x);

++N_;
}
	  }	// namespace jtutil::

//******************************************************************************

//
// these functions compute the corresponding norms
//
namespace jtutil
	  {
template<typename fp_t>
  fp_t norm<fp_t>::mean() const { return sum_/fp_t(N_); }
template<typename fp_t>
  fp_t norm<fp_t>::std_dev() const
        {
        if (is_empty()) return fp_t(0);
        return sqrt(jtutil::max(fp_t(0),fp_t(N_)*sum2_-sum_*sum_))/fp_t(N_);
        }
template<typename fp_t>
  fp_t norm<fp_t>::two_norm() const { return sqrt(sum2_); }
template<typename fp_t>
  fp_t norm<fp_t>::rms_norm() const
	{ assert(is_nonempty()); return sqrt(sum2_/fp_t(N_)); }
	  }	// namespace jtutil::

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// ***** template instantiations *****
//

template class jtutil::norm<float>;
template class jtutil::norm<double>;
