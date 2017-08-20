// cpm_map.cc -- "integer +/-" mapping  i --> j = const +/- i
// $Header$
//
// jtutil::cpm_map::cpm_map  # mirror map, specified by fixed point
// jtutil::cpm_map::cpm_map  # generic map specified by sample point & sign
// jtutil::cpm_map::cpm_map  # generic map specified by *fp* sample point & sign
//
// ***** template instantiations *****
//

#include <assert.h>
#include <stdio.h>
#include <cmath>

#include "../cctk.h"
#include "util.hh"
#include "cpm_map.hh"

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function constructs a  cpm_map  object with a "-" sign and a
// specified fixed point (must be integer or half-integer) and domain.
// The sample point need not be in the map's domain/range.
//
namespace jtutil
	  {
template <typename fp_t>
cpm_map<fp_t>::cpm_map(int min_i_in, int max_i_in,
		       fp_t fixed_point)
	: min_i_(min_i_in), max_i_(max_i_in),
	  map_is_plus_(false)
{
const fp_t d_offset = 2.0 * fixed_point;
if (! fuzzy<fp_t>::is_integer(d_offset))
   then error_exit(ERROR_EXIT,
"***** cpm_map::cpm_map (mirror):\n"
"        fixed_point=%g isn't (fuzzily) integral or half-integral!\n"
,
		   double(fixed_point));			/*NOTREACHED*/

offset_ = round<fp_t>::to_integer(d_offset);

// verify that we have setup correct
assert(
   map_unchecked(fuzzy<fp_t>::floor  (fixed_point))
   ==
		 fuzzy<fp_t>::ceiling(fixed_point)
      );
}
	  }	// namespace jtutil::

//******************************************************************************

//
// This function constructs a generic  cpm_map  object, with the mapping
// specified by a sample point  sample_i --> sample_j  and by sign.
// The sample point need not be in the map's domain/range.
//
namespace jtutil
	  {
template <typename fp_t>
cpm_map<fp_t>::cpm_map(int min_i_in, int max_i_in,
		       int sample_i, int sample_j,
		       bool map_is_plus_in)
	: min_i_(min_i_in), max_i_(max_i_in),
	  offset_(map_is_plus_in ? sample_j - sample_i
			  : sample_j + sample_i),
	  map_is_plus_(map_is_plus_in)
{
// verify that we have setup correct
assert( map_unchecked(sample_i) == sample_j );
}
	  }	// namespace jtutil::

//******************************************************************************

//
// This function constructs a generic  cpm_map  object, with the mapping
// specified by a *fp* sample point  sample_i --> sample_j  (which
// must specify an  integer --> integer  mapping, i.e.  4.2 --> 4.2  is
// ok for a + map, and 4.5 --> 4.5 is ok for a minus map, but  4.2 --> 4.7
// is never ok) and by sign.  The sample point need not be in the map's
// domain/range.
//
namespace jtutil
	  {
template <typename fp_t>
cpm_map<fp_t>::cpm_map(int min_i_in, int max_i_in,
		       fp_t sample_i, fp_t sample_j,
		       bool map_is_plus_in)
	: min_i_(min_i_in), max_i_(max_i_in),
	  map_is_plus_(map_is_plus_in)
{
const fp_t fp_offset = map_is_plus_in ? sample_j - sample_i
				    : sample_j + sample_i;
if (! fuzzy<fp_t>::is_integer(fp_offset))
   then error_exit(ERROR_EXIT,
"***** cpm_map::cpm_map (generic via fp sample point):\n"
"        fp_offset=%g isn't fuzzily integral!\n"
"        ==> sample_i=%g --> sample_j=%g\n"
"            doesn't fuzzily specify an  integer --> integer  mapping!\n"
,
		   double(fp_offset),
		   double(sample_i), double(sample_j));		/*NOTREACHED*/

offset_ = round<fp_t>::to_integer(fp_offset);

// verify that we have setup correct
assert(
   map_unchecked(      fuzzy<fp_t>::floor(sample_i))
   ==
   (map_is_plus_in ? fuzzy<fp_t>::floor  (sample_j)
		   : fuzzy<fp_t>::ceiling(sample_j)) );
}
	  }	// namespace jtutil::

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// ***** template instantiations *****
//

template class jtutil::cpm_map<float>;
template class jtutil::cpm_map<double>;
