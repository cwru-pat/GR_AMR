#ifndef AHFD_ARRAY_HH
#define AHFD_ARRAY_HH

//
// prerequisites:
//    <assert.h>
//    <stddef.h>	// for NULL
//    "stdc.h"
//    "util.hh"		// for jtutil::how_many_in_range()
//

//
// The templates defined in this file represent n-dimensional row-major
// contiguous arrays, parameterized by the data type (most commonly float
// or double, but could also be bool, int, long double, ...).  The
// underlying storage can be either supplied by the caller, or allocated
// by new[].  In the latter case, arbitrary strides are also possible.
// Unfortunately, allowing arbitrary strides makes subscripting a bit
// more expensive. :(
//
// These arrays cannot be copied or passed to functions by value; use
// pass by reference instead.  (This is a feature, not a bug: passing
// large arrays by value would be ++slow, so we don't want to run the
// risk of this happening accidentally.)
//

//
// Stroustrup ("The C++ Programming Language", 3rd edition, appendix C.7)
// suggests the use of STL vectors of STL vectors to provide multidimensional
// arrays.  However, those "arrays" aren't contiguous in memory, so the
// compiler can't do strength reduction on any but the last subscript
// when the arrays are accessed in a loop.  In contrast, the arrays
// defined here are contiguous, and all subscripts can (should, if the
// compiler is good) be strength-reduced.
//
// The STL valarray templates offer a superset of the functionality of
// these templates, but they've only recently been provided with gcc;
// at some time in the future I may migrate these classes to become
// valarray wrappers.
//
// Or, even better, boost (http://www.boost.org) is working on a
// multi_array<> template, which looks very nice, abeit complicated.
// See
//   http://groups.yahoo.com/group/boost/message/10230
// for the development discussions for this.
//

//******************************************************************************

namespace jtutil
	  {
template <typename T>
class	array1d
	{
public:
	// array info
	int min_i() const { return min_i_; }
	int max_i() const { return max_i_; }
	int N_i() const { return jtutil::how_many_in_range(min_i_, max_i_); }
	bool is_valid_i(int i) const { return (i >= min_i_) && (i <= max_i_); }

	// subscripting functions
	// (low-level, dangerous, use with caution!)
	// FIXME: should we also provide the reverse mapping, i.e.
	//        subscript --> (i) ?
	int subscript_unchecked(int i) const
		{ return offset_ + stride_i_*i; }
	int subscript(int i) const
		{
		// n.b. we want each assert() here to be on a separate
		//	source line, so an assert() failure message can
		//	pinpoint *which* index is bad
		assert( is_valid_i(i) );
		const int posn = subscript_unchecked(i);
		assert(posn >= 0);
		assert(posn <= max_subscript_);
		return posn;
		}
	int subscript_offset() const { return offset_; }
	int subscript_stride_i() const { return stride_i_; }
		

	// normal-use access functions
	// ... rvalue
	const T& operator()(int i) const { return array_[ subscript(i) ]; }
	// ... lvalue
	      T& operator()(int i)       { return array_[ subscript(i) ]; }

	// get access to internal 0-origin 1D storage array
	// (low-level, dangerous, use with caution!)
	// ... semantics of N_array() may not be what you want
	//     if strides specify noncontiguous storage
	int N_array() const { return max_subscript_+stride_i_; }
	const T* data_array() const { return const_cast<const T*>(array_); }
	      T* data_array()       { return array_; }

	// constructor, destructor
	// ... constructor initializes all array elements to T(0.0)
	// ... omitted strides default to C storage order
	array1d(int min_i_in, int max_i_in,
		T *array_in = NULL,	// caller-provided storage array
					// if non-NULL
		int stride_i_in = 0);
	~array1d();

private:
	// we forbid copying and passing by value
	// by declaring the copy constructor and assignment operator
	// private, but never defining them
	array1d(const array1d<T>& rhs);
	array1d<T>& operator=(const array1d<T>& rhs);

private:
	// n.b. we declare the array pointer first in the class
	// ==> it's probably at 0 offset
	// ==> we may get slightly faster array access
	T* array_;		// --> new-allocated 1D storage array

	// subscripting info
	// n.b. put this next in class so it should be in the same
	//	cpu cache line as  array_  ==> faster array access
	int offset_, stride_i_;

	// min/max array bounds
	const int min_i_, max_i_;
	int max_subscript_;

	// n.b. put this at end of class since performance doesn't matter
	bool we_own_array_;	// true ==> array_ --> new[] array which we own
				// false ==> array_ --> client-owned storage
	};

//******************************************************************************


template <typename T>
class	array2d
	{
public:
	// array info
	int min_i() const { return min_i_; }
	int max_i() const { return max_i_; }
	int min_j() const { return min_j_; }
	int max_j() const { return max_j_; }
	int N_i() const { return jtutil::how_many_in_range(min_i_, max_i_); }
	int N_j() const { return jtutil::how_many_in_range(min_j_, max_j_); }
	bool is_valid_i(int i) const { return (i >= min_i_) && (i <= max_i_); }
	bool is_valid_j(int j) const { return (j >= min_j_) && (j <= max_j_); }
	bool is_valid_ij(int i, int j) const
		{ return is_valid_i(i) && is_valid_j(j); }

	// subscripting functions
	// (low-level, dangerous, use with caution!)
	// FIXME: should we also provide the reverse mapping, i.e.
	//        subscript --> (i,j) ?
	int subscript_unchecked(int i, int j) const
		{ return offset_ + stride_i_*i + stride_j_*j; }
	int subscript(int i, int j) const
		{
		// n.b. we want each assert() here to be on a separate
		//	source line, so an assert() failure message can
		//	pinpoint *which* index is bad
		assert( is_valid_i(i) );
		assert( is_valid_j(j) );
		const int posn = subscript_unchecked(i,j);
		assert(posn >= 0);
		assert(posn <= max_subscript_);
		return posn;
		}
	int subscript_offset() const { return offset_; }
	int subscript_stride_i() const { return stride_i_; }
	int subscript_stride_j() const { return stride_j_; }

	// normal-use access functions
	// ... rvalue
	const T& operator()(int i, int j) const
		{ return array_[ subscript(i,j) ]; }
	// ... lvalue
	      T& operator()(int i, int j)
		{ return array_[ subscript(i,j) ]; }

	// get access to internal 0-origin 1D storage array
	// (low-level, dangerous, use with caution!)
	// ... semantics of N_array() may not be what you want
	//     if strides specify noncontiguous storage
	int N_array() const { return max_subscript_+stride_j_; }
	const T* data_array() const { return const_cast<const T*>(array_); }
	      T* data_array()       { return array_; }

	// constructor, destructor
	// ... constructor initializes all array elements to T(0.0)
	// ... omitted strides default to C storage order
	array2d(int min_i_in, int max_i_in,
		int min_j_in, int max_j_in,
		T *array_in = NULL,	// caller-provided storage array
					// if non-NULL
		int stride_i_in = 0, int stride_j_in = 0);
	~array2d();

private:
	// we forbid copying and passing by value
	// by declaring the copy constructor and assignment operator
	// private, but never defining them
	array2d(const array2d<T>& rhs);
	array2d<T>& operator=(const array2d<T>& rhs);

private:
	// n.b. we declare the array pointer first in the class
	// ==> it's probably at 0 offset
	// ==> we may get slightly faster array access
	T* array_;		// --> new-allocated 1D storage array

	// subscripting info
	// n.b. put this next in class so it should be in the same
	//	cpu cache line as  array_  ==> faster array access
	int offset_, stride_i_, stride_j_;

	// min/max array bounds
	const int min_i_, max_i_;
	const int min_j_, max_j_;
	int max_subscript_;

	// n.b. put this at end of class since performance doesn't matter
	bool we_own_array_;	// true ==> array_ --> new[] array which we own
				// false ==> array_ --> client-owned storage
	};

//******************************************************************************


template <typename T>
class	array3d
	{
public:
	// array info
	int min_i() const { return min_i_; }
	int max_i() const { return max_i_; }
	int min_j() const { return min_j_; }
	int max_j() const { return max_j_; }
	int min_k() const { return min_k_; }
	int max_k() const { return max_k_; }
	int N_i() const { return jtutil::how_many_in_range(min_i_, max_i_); }
	int N_j() const { return jtutil::how_many_in_range(min_j_, max_j_); }
	int N_k() const { return jtutil::how_many_in_range(min_k_, max_k_); }
	bool is_valid_i(int i) const { return (i >= min_i_) && (i <= max_i_); }
	bool is_valid_j(int j) const { return (j >= min_j_) && (j <= max_j_); }
	bool is_valid_k(int k) const { return (k >= min_k_) && (k <= max_k_); }
	bool is_valid_ijk(int i, int j, int k) const
		{ return is_valid_i(i) && is_valid_j(j) && is_valid_k(k); }

	// subscripting functions
	// (low-level, dangerous, use with caution!)
	// FIXME: should we also provide the reverse mapping, i.e.
	//        subscript --> (i,j,k) ?
	int subscript_unchecked(int i, int j, int k) const
		{ return offset_ + stride_i_*i + stride_j_*j + stride_k_*k; }
	int subscript(int i, int j, int k) const
		{
		// n.b. we want each assert() here to be on a separate
		//	source line, so an assert() failure message can
		//	pinpoint *which* index is bad
		assert( is_valid_i(i) );
		assert( is_valid_j(j) );
		assert( is_valid_k(k) );
		const int posn = subscript_unchecked(i,j,k);
		assert(posn >= 0);
		assert(posn <= max_subscript_);
		return posn;
		}
	int subscript_offset() const { return offset_; }
	int subscript_stride_i() const { return stride_i_; }
	int subscript_stride_j() const { return stride_j_; }
	int subscript_stride_k() const { return stride_k_; }

	// normal-use access functions
	// ... rvalue
	const T& operator()(int i, int j, int k) const
		{ return array_[ subscript(i,j,k) ]; }
	// ... lvalue
	      T& operator()(int i, int j, int k)
		{ return array_[ subscript(i,j,k) ]; }

	// get access to internal 0-origin 1D storage array
	// (low-level, dangerous, use with caution!)
	// ... semantics of N_array() may not be what you want
	//     if strides specify noncontiguous storage
	int N_array() const { return max_subscript_+stride_k_; }
	const T* data_array() const { return const_cast<const T*>(array_); }
	      T* data_array()       { return array_; }

	// constructor, destructor
	// ... constructor initializes all array elements to T(0.0)
	// ... omitted strides default to C storage order
	array3d(int min_i_in, int max_i_in,
		int min_j_in, int max_j_in,
		int min_k_in, int max_k_in,
		T *array_in = NULL,	// caller-provided storage array
					// if non-NULL
		int stride_i_in = 0, int stride_j_in = 0, int stride_k_in = 0);
	~array3d();

private:
	// we forbid copying and passing by value
	// by declaring the copy constructor and assignment operator
	// private, but never defining them
	array3d(const array3d<T>& rhs);
	array3d<T>& operator=(const array3d<T>& rhs);

private:
	// n.b. we declare the array pointer first in the class
	// ==> it's probably at 0 offset
	// ==> we may get slightly faster array access
	T* array_;		// --> new-allocated 1D storage array

	// subscripting info
	// n.b. put this next in class so it should be in the same
	//	cpu cache line as  array_  ==> faster array access
	int offset_, stride_i_, stride_j_, stride_k_;

	// min/max array bounds
	const int min_i_, max_i_;
	const int min_j_, max_j_;
	const int min_k_, max_k_;
	int max_subscript_;

	// n.b. put this at end of class since performance doesn't matter
	bool we_own_array_;	// true ==> array_ --> new[] array which we own
				// false ==> array_ --> client-owned storage
	};

//******************************************************************************

	  }	// namespace jtutil::

//******************************************************************************

#endif	/* AHFINDERDIRECT__ARRAY_HH */
