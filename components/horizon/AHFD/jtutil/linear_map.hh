// linear_map.hh -- linear mapping from integers <--> floating point values
// $Header$
//
// jtutil::linear_map - linear mapping from integers <--> floating point values
//

#ifndef AHFINDERDIRECT__LINEAR_MAP_HH
#define AHFINDERDIRECT__LINEAR_MAP_HH

//
// prerequisites
//    <assert.h>
//    "stdc.h"
//    "util.hh"		// for jtutil::how_many_in_range() and fuzzy::
//

//
// The template defined in this file represents a linear mapping between
// a contiguous range of integers, and an arithmetic progression of
// floating point values, parameterized by the floating point data type.
//

#ifndef NDEBUG
//
// Full bounds checking is done on the mapping in both directions
//
#endif

//*****************************************************************************

namespace jtutil
	  {
template <typename fp_t>
class	linear_map
	{
public:
	// integer bounds info
	int min_int() const { return min_int_; }
	int max_int() const { return max_int_; }
	int N_points() const
		{ return jtutil::how_many_in_range(min_int_,max_int_); }
	bool is_in_range(int i) const
		{ return (i >= min_int()) && (i <= max_int()); }
	int clamp(int i) const
		{
		if	(i < min_int())
		   then return min_int();
		else if (i > max_int())
		   then return max_int();
		else	return i;
		}

	// convert int --> fp
	fp_t fp_of_int_unchecked(int i) const
		{ return origin_ + delta_*i; }
	fp_t fp_of_int(int i) const
		{
		assert(is_in_range(i));
		return fp_of_int_unchecked(i);
		}

	// converg delta_int --> delta_fp
	fp_t delta_fp_of_delta_int(int delta_i) const
		{ return delta_ * delta_i; }

	// fp bounds info
	fp_t origin() const { return origin_; }
	fp_t delta_fp() const { return delta_; }
	fp_t inverse_delta_fp() const { return inverse_delta_; }
	fp_t min_fp() const { return fp_of_int_unchecked(min_int_); }
	fp_t max_fp() const { return fp_of_int_unchecked(max_int_); }
	bool is_in_range(fp_t x) const
		{
		return    fuzzy<fp_t>::GE(x,min_fp())
		       && fuzzy<fp_t>::LE(x,max_fp());
		}
	fp_t clamp(fp_t x) const
		{
		if	(x < min_fp())
		   then return min_fp();
		else if (x > max_fp())
		   then return max_fp();
		else	return x;
		}

	// convert linear map indices <--> C-style 0-origin indices
	int zero_origin_int(int i) const { return i - min_int(); }
	int map_int(int zero_origin_i) { return zero_origin_i + min_int(); }

	// convert fp --> int coordinate, but return result as fp
	// (which need not be fuzzily integral)
	fp_t fp_int_of_fp(fp_t x) const;

	// convert fp --> int, check being fuzzily integral
	enum	noninteger_action	// what to do if "int"
					// isn't fuzzily integral?
		{
		nia_error,		// jtutil::error_exit(...)
		nia_warning,		// print warning msg,
					// then round to nearest
		nia_round,		// (silently) round to nearest
		nia_floor,		// (silently) round to -infinity
		nia_ceiling		// (silently) round to +infinity
		};
	int int_of_fp(fp_t x, noninteger_action nia = nia_error) const;

	// convert delta_fp --> delta_int, check being fuzzily integral
	int delta_int_of_delta_fp(fp_t delta_x,
				  noninteger_action nia = nia_error)
		const;

	// constructors
	linear_map(int min_int_in, int max_int_in,
		   fp_t min_fp_in, fp_t delta_fp_in, fp_t max_fp_in);
	// ... construct with subrange of existing linear_map
	linear_map(const linear_map<fp_t> &lm_in,
		   int min_int_in, int max_int_in);

	// no need for explicit destructor, compiler-generated no-op is ok

	// no need for copy constructor or assignment operator,
	// compiler-generated defaults are ok

private:
	// common code (argument validation & setup) for all constructors
	// assumes min_int_, max_int_, delta_ already initialized,
	//         other class members *not* initialized
	void constructor_common(fp_t min_fp_in, fp_t max_fp_in);

	// these define the actual mapping
	// via the  fp_of_int()  function (above)
	fp_t origin_, delta_;

	// cache of 1.0/delta_
	// ==> avoids fp divide in inverse_delta_fp()
	// ==> also makes fp --> int conversions slightly faster
	fp_t inverse_delta_;

	// bounds (inclusive)
	const int min_int_, max_int_;
	};
	  }	// namespace jtutil::

//******************************************************************************

#endif	/* AHFINDERDIRECT__LINEAR_MAP_HH */
