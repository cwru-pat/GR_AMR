// linear_map.cc -- linear mapping from integers <--> floating point values
// $Header$
//
// jtutil::linear_map::linear_map - basic ctor
// jtutil::linear_map::linear_map - ctor for subrange of existing linear_map
// jtutil::linear_map::fp_int_of_fp - convert fp --> int coord, return as fp
// jtutil::linear_map::int_of_fp - convert fp --> int, check for ==fuzzy int
// jtutil::linear_map::delta_int_of_delta_fp - ... same for spacings
//
// ***** template instantiation *****
//

#include <assert.h>
#include <stdio.h>
#include <cmath>

#include "../cctk.h"
#include "util.hh"
#include "linear_map.hh"

//******************************************************************************
//******************************************************************************
// linear_map -- linear mapping from integers <--> floating point values
//******************************************************************************
//******************************************************************************

//
// This function constructs a  linear_map<fp_t>  object.
//
namespace jtutil
	  {
template <typename fp_t>
linear_map<fp_t>::linear_map(int min_int_in, int max_int_in,
			     fp_t min_fp_in, fp_t delta_fp_in, fp_t max_fp_in)
	: delta_(delta_fp_in), inverse_delta_(1.0 / delta_fp_in),
	  min_int_(min_int_in), max_int_(max_int_in)
{
constructor_common(min_fp_in, max_fp_in);
}
	  }	// namespace jtutil::

//******************************************************************************

//
// This function constructs a  linear_map<fp_t>  object with a subrange
// of an existing one.
//
namespace jtutil
	  {
template <typename fp_t>
linear_map<fp_t>::linear_map(const linear_map<fp_t> &lm_in,
			     int min_int_in, int max_int_in)	// subrange
	: delta_(lm_in.delta_fp()), inverse_delta_(lm_in.inverse_delta_fp()),
	  min_int_(min_int_in), max_int_(max_int_in)
{
if (! (is_in_range(min_int_in) && is_in_range(max_int_in)) )
   then error_exit(ERROR_EXIT,
"***** linear_map<fp_t>::linear_map:\n"
"        min_int_in=%d and/or max_int_in=%d\n"
"        aren't in integer range [%d,%d] of existing linear_map!\n"
,
		   min_int_, max_int_,
		   lm_in.min_int(), lm_in.max_int());		/*NOTREACHED*/

constructor_common(lm_in.fp_of_int_unchecked(min_int_in),
		   lm_in.fp_of_int_unchecked(max_int_in));
}
	  }	// namespace jtutil::

//******************************************************************************

//
// This function does the common argument validation and setup for
// all the constructors of class  linear_map<fp_t>:: .
//
namespace jtutil
	  {
template <typename fp_t>
void linear_map<fp_t>::constructor_common(fp_t min_fp_in, fp_t max_fp_in)
	// assumes
	//	min_int_, max_int_, delta_, inverse_delta_
	// are already initialized
	// ==> ok to use min_int(), max_int(), delta_fp(), inverse_delta_fp()
	// ... other class members *not* yet initialized
{
origin_ = 0.0;		// temp value
origin_ = min_fp_in - fp_of_int_unchecked(min_int());

// this should be guaranteed by the above calculation
assert( fuzzy<fp_t>::EQ(fp_of_int_unchecked(min_int()), min_fp_in) );

// this is a test of the consistency of the input arguments
if (fuzzy<fp_t>::NE(fp_of_int_unchecked(max_int()), max_fp_in))
   then error_exit(ERROR_EXIT,
"***** linear_map<fp_t>::linear_map:\n"
"        int range [%d,%d]\n"
"        and fp range [%g(%g)%g]\n"
"        are (fuzzily) inconsistent!\n"
,
		   min_int(), max_int(),
		   double(min_fp_in), double(delta_fp()), double(max_fp_in));
								/*NOTREACHED*/
}
	  }	// namespace jtutil::

//******************************************************************************

//
// This function converts  fp  --> int  coordinate, returning the result
// as an fp (which need not be fuzzily integral).
//
namespace jtutil
	  {
template <typename fp_t>
fp_t linear_map<fp_t>::fp_int_of_fp(fp_t x)
	const
{
if (! is_in_range(x))
   then error_exit(ERROR_EXIT,
"***** linear_map<fp_t>::fp_int_of_fp:\n"
"        fp value x=%g is (fuzzily) outside the grid!\n"
"        {min(delta)max}_fp = %g(%g)%g\n"
,
		   double(x),
		   double(min_fp()), double(delta_fp()), double(max_fp()));
								/*NOTREACHED*/

return inverse_delta_ * (x - origin_);
}
	  }	// namespace jtutil::

//******************************************************************************

//
// This function converts  fp  --> int  and checks that the result is
// fuzzily integral.  (The  nia  argument specifies what to do if the
// result *isn't* fuzzily integral.)
//
// FIXME:
// Having to explicitly specify the namespace for jtutil::round<fp_t>::
// is ++ugly. :(
//
namespace jtutil
	  {
template <typename fp_t>
int linear_map<fp_t>::int_of_fp(fp_t x, noninteger_action nia /* = nia_error */)
	const
{
const fp_t fp_int = fp_int_of_fp(x);

if (fuzzy<fp_t>::is_integer(fp_int))
   then {
	// x is (fuzzily) a grid point ==> return that
	return jtutil::round<fp_t>::to_integer(fp_int);	// *** EARLY RETURN ***
	}

// get to here ==> x isn't (fuzzily) a grid point
static const char *const noninteger_msg =
   "%s linear_map<fp_t>::int_of_fp:\n"
   "        x=%g isn't (fuzzily) a grid point!\n"
   "        {min(delta)max}_fp() = %g(%g)%g\n"
   ;
switch	(nia)
	{
case nia_error:
	error_exit(ERROR_EXIT,
		   noninteger_msg,
		   "*****",
		   double(x),
		   double(min_fp()), double(delta_fp()), double(max_fp()));
								/*NOTREACHED*/

case nia_warning:
	printf(noninteger_msg,
	       "---",
	       double(x),
	       double(min_fp()), double(delta_fp()), double(max_fp()));
	// fall through

case nia_round:
	return jtutil::round<fp_t>::to_integer(fp_int);	// *** EARLY RETURN ***

case nia_floor:
	return jtutil::round<fp_t>::floor(fp_int);	// *** EARLY RETURN ***

case nia_ceiling:
	return jtutil::round<fp_t>::ceiling(fp_int);	// *** EARLY RETURN ***

default:
	error_exit(PANIC_EXIT,
"***** linear_map<fp_t>::int_of_fp: illegal nia=(int)%d\n"
"                                   (this should never happen!)\n"
,
		   int(nia));					/*NOTREACHED*/
	}
return 0;	// dummy return to quiet gcc
		// (which doesn't grok that error_exit() never returns)
}
	  }	// namespace jtutil::

//******************************************************************************

//
// This function converts "delta" spacings in the fp coordinate to
// corresponding "delta" spacings in the int coordinate, and checks that
// the result is fuzzily integral.  (The  nia  argument specifies what to
// do if the result *isn't* fuzzily integral.)
//
// FIXME:
// Having to explicitly specify the namespace for jtutil::round<fp_t>::
// is ++ugly. :(
//
namespace jtutil
	  {
template <typename fp_t>
int linear_map<fp_t>::delta_int_of_delta_fp
	(fp_t delta_x, noninteger_action nia /* = nia_error */)
	const
{
const fp_t fp_delta_int = inverse_delta_ * delta_x;

if (fuzzy<fp_t>::is_integer(fp_delta_int))
   then {
	// delta_x is (fuzzily) an integer number of grid spacings
	// ==> return that
	return jtutil::round<fp_t>::to_integer(fp_delta_int);
							// *** EARLY RETURN ***
	}

// get to here ==> delta_x isn't (fuzzily) an integer number of grid spacings
static const char *const noninteger_msg =
   "%s linear_map<fp_t>::delta_int_of_delta_fp:\n"
   "        delta_x=%g isn't (fuzzily) an integer number of grid spacings!\n"
   "        {min(delta)max}_fp() = %g(%g)%g\n"
   ;
switch	(nia)
	{
case nia_error:
	error_exit(ERROR_EXIT,
		   noninteger_msg,
		   "*****",
		   double(delta_x),
		   double(min_fp()), double(delta_fp()), double(max_fp()));
								/*NOTREACHED*/

case nia_warning:
	printf(noninteger_msg,
	       "---",
	       double(delta_x),
	       double(min_fp()), double(delta_fp()), double(max_fp()));
	// fall through

case nia_round:
	return jtutil::round<fp_t>::to_integer(fp_delta_int);
							// *** EARLY RETURN ***

case nia_floor:
	return jtutil::round<fp_t>::floor(fp_delta_int);// *** EARLY RETURN ***

case nia_ceiling:
	return jtutil::round<fp_t>::ceiling(fp_delta_int);
							// *** EARLY RETURN ***

default:
	error_exit(PANIC_EXIT,
"***** linear_map<fp_t>::delta_int_of_delta_fp: illegal nia=(int)%d\n"
"                                               (this should never happen!)\n"
,
		   int(nia));					/*NOTREACHED*/
	}
return 0;	// dummy return to quiet gcc
		// (which doesn't grok that error_exit() never returns)
}
	  }	// namespace jtutil::

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// ***** template instantiation *****
//

template class jtutil::linear_map<float>;
template class jtutil::linear_map<double>;
