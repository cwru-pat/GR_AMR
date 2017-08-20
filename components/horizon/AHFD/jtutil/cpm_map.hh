// cpm_map.hh -- "integer +/-" mapping  i --> j = const +/- i
// $Header$
//
// cpm_map - "integer +/-" mapping
//

#ifndef AHFINDERDIRECT__CPM_MAP_HH
#define AHFINDERDIRECT__CPM_MAP_HH

//
// prerequisites
//	<assert.h>
//	"stdc.h"
//	"util.hh"		// for jtutil::how_many_in_range()
//

//
// This class represents a unit-slope linear integer mapping
//	i --> j = const +/- i
//

#ifndef NDEBUG
//
// Full bounds checking is done on the mapping in both directions.
//
#endif

//******************************************************************************

namespace jtutil
	  {
template <typename fp_t>
class	cpm_map
	{
public:
	// bounds info -- domain
	int min_i() const { return min_i_; }
	int max_i() const { return max_i_; }
	int N_points() const
		{ return jtutil::how_many_in_range(min_i_,max_i_); }
	bool in_domain(int i) const { return (i >= min_i_) && (i <= max_i_); }

	// is the mapping + or - ?
	bool is_plus() const { return map_is_plus_; }
	bool is_minus() const { return !map_is_plus_; }
	int sign() const { return map_is_plus_ ? +1 : -1; }
	fp_t fp_sign() const { return map_is_plus_ ? +1.0 : -1.0; }

	// the mapping itself
	int map_unchecked(int i) const
		{
		return map_is_plus_ ? offset_ + i
				    : offset_ - i;
		}
	int inv_map_unchecked(int j) const
		{
		return map_is_plus_ ? j - offset_
				    : offset_ - j;
		}
	int map(int i) const
		{
		assert(in_domain(i));
		return map_unchecked(i);
		}
	int inv_map(int j) const
		{
		int i = inv_map_unchecked(j);
		assert(in_domain(i));
		return i;
		}

	// bounds info -- range
	// ... we use the unchecked map here in case the domain is empty
	int min_j() const
		{
		return map_is_plus_ ? map_unchecked(min_i_)
				    : map_unchecked(max_i_);
		}
	int max_j() const
		{
		return map_is_plus_ ? map_unchecked(max_i_)
				    : map_unchecked(min_i_);
		}
	bool in_range(int j) const { return in_domain(inv_map_unchecked(j)); }

	//
	// constructors
	//

	// "mirror" map: i --> const - i
	// ... map specified by fixed point (must be integer or half-integer)
	// ... fixed point need not be in domain/range
	cpm_map(int min_i_in, int max_i_in,
		fp_t fixed_point);

	// "shift" map: i --> const + i
	// ... map specified by shift amount
	// ... default is identity map
	cpm_map(int min_i_in, int max_i_in,
		int shift_amount = 0)
		: min_i_(min_i_in), max_i_(max_i_in),
		  offset_(shift_amount), map_is_plus_(true)
		{ }

	// generic map: i --> const +/- i
	// ... map specified by sample point sample_i --> sample_j
	//     and by sign (one of  {plus,minus}_map )
	// ... sample point need not be in domain/range
	cpm_map(int min_i_in, int max_i_in,
		int sample_i, int sample_j,
		bool map_is_plus_in);

	// generic map: i --> const +/- i
	// ... map specified by *fp* sample point sample_i --> sample_j
	//     (must specify an integer --> integer mapping)
	//     and by sign (one of  {plus,minus}_map )
	// ... hence if sign is -1, then sample_i and sample_j
	//     must both be half-integral
	// ... sample point need *not* be in domain/range
	cpm_map(int min_i_in, int max_i_in,
		fp_t sample_i, fp_t sample_j,
		bool map_is_plus_in);

	// no need for explicit destructor, compiler-generated no-op is ok
	// ditto for copy constructor and assignment operator

private:
	// bounds (inclusive)
	int min_i_, max_i_;

	// these define the actual mapping
	int offset_;
	bool map_is_plus_;
	};
	  }	// namespace jtutil::

//******************************************************************************

#endif	/* AHFINDERDIRECT__CPM_MAP_HH */
