#ifndef AHFD_UTIL_HH
#define AHFD_UTIL_HH

#include <string>
#include "../../../../cosmo_macros.h"
namespace jtutil
	  {

//******************************************************************************

// how many integers are in the closed interval [low,high]
inline int how_many_in_range(int low, int high) { return high - low + 1; }

// is an integer even/odd
inline int is_even(int i) { return !(i & 0x1); }
inline int is_odd (int i) { return  (i & 0x1); }

//
// min/max/absolute value template
// FIXME: <algorithm> is supposed to have min/max, but these are
//	  broken on too many platforms (eg older gcc versions)
//
template <typename T>
  inline T min(T x, T y) { return (x < y) ? x : y; }
template <typename T>
  inline T max(T x, T y) { return (x > y) ? x : y; }
template <typename T>
  inline T abs(T x) { return (x > 0) ? x : -x; }

//
// These functions raise their arguments to various small-integer powers.
//
template <typename T>
  inline T pow2(T x) { return x*x; }
template <typename T>
  inline T pow3(T x) { return x*x*x; }
template <typename T>
  inline T pow4(T x) { return pow2(pow2(x)); }
#ifdef NOT_USED
template <typename T>
  inline T pow5(T x) { return x * pow4(x); }
template <typename T>
  inline T pow6(T x) { return pow3(pow2(x)); }
template <typename T>
  inline T pow7(T x) { return x * pow6(x); }
template <typename T>
  inline T pow8(T x) { return pow2(pow2(pow2(x))); }
#endif

//
// misc math stuff
//
double signum(double x);
double hypot3(double x, double y, double z);
double arctan_xy(double x, double y);

// reduce x modulo xmod to be in the interval [xmin,xmax],
// or error_exit() if no such value exists
double modulo_reduce(double x, double xmod, double xmin, double xmax);

// initialize C-style array to all zeros
template <typename fp_t>
  void zero_C_array(int N, fp_t array[]);

//
// more misc math stuff, valid only if <math.h> has been #included;
//
#ifdef M_PI	/* test for <math.h> */
		/* n.b. C-style comment needed for some preprocessors! */
  // floor/ceiling of double, returned as an int
  inline int ifloor(double x) { return static_cast<int>(floor(x)); }
  inline int iceil (double x) { return static_cast<int>(ceil (x)); }
#endif

#ifdef PI	/* PI is defined in "../include/stdc.h" */
		/* n.b. C-style comment needed for some preprocessors! */
  // convert degrees <--> radians
  template <typename fp_t>
    inline fp_t degrees_of_radians(fp_t radians) { return (180.0/PI)*radians; }
  template <typename fp_t>
    inline fp_t radians_of_degrees(fp_t degrees) { return (PI/180.0)*degrees; }
#endif

//******************************************************************************

//
// This template class computes means, 2-norms, rms-norms, infinity-norms,
// and various related values.
//
template <typename fp_t>
class	norm
	{
public:
	// get norms etc
	fp_t mean() const;
	fp_t std_dev() const;         // sqrt(sum (x_i - average of x_i)^2 / N)
	fp_t two_norm() const;		// sqrt(sum x_i^2)
	fp_t rms_norm() const;		// sqrt(average of x_i^2)
	fp_t infinity_norm() const { return max_abs_value_; }

	fp_t max_abs_value() const { return max_abs_value_; }
	fp_t min_abs_value() const { return min_abs_value_; }

	fp_t max_value() const { return max_value_; }
	fp_t min_value() const { return min_value_; }

	// specify data point
	void data(fp_t x);

	// have any data points been specified?
	bool is_empty()    const { return N_ == 0; }
	bool is_nonempty() const { return N_ > 0; }

	// reset ==> just like newly-constructed object
	void reset();

	// constructor, destructor
	// ... compiler-generated no-op destructor is ok
	norm();

private:
        // we forbid copying and passing by value
        // by declaring the copy constructor and assignment operator
        // private, but never defining them
        norm(const norm &rhs);
	norm& operator=(const norm &rhs);

private:
	long N_;		// # of data points
	fp_t sum_;		// sum(data)
	fp_t sum2_;		// sum(data^2)
	fp_t max_abs_value_;	// max |data|
	fp_t min_abs_value_;	// min |data|
	fp_t max_value_;	// max data
	fp_t min_value_;	// min data
	};

//******************************************************************************

//
// This template does fuzzy comparisons and related operations on
// floating point values, parameterized by the floating point type.
//
// The fuzzy comparison semantics are based on those of APL, but are
// modified to use an absolute error tolerance for values close to 0.
//

// this template class has only static members
// ... it's a class, not a namespace, because we want to express the
//     semantics that the entire class is a single template, rather
//     than the individual members being conceptually-unrelated templates
// ... moreover, we need the *data* member (template)  tolerance , and
//     it seems C++ doesn't grok data templates which aren't in a class
template <typename fp_t>
class	fuzzy
	{
public:
	// comparison tolerance (may be modified by user code if needed)
	static fp_t get_tolerance() { return tolerance_; }
	static void set_tolerance(fp_t new_tolerance)
		{ tolerance_ = new_tolerance; }

	// fuzzy commparisons
	static bool EQ(fp_t x, fp_t y);
	static bool NE(fp_t x, fp_t y) { return ! EQ(x,y); }
	static bool LT(fp_t x, fp_t y) { return EQ(x,y) ? false : (x < y); }
	static bool LE(fp_t x, fp_t y) { return EQ(x,y) ? true  : (x < y); }
	static bool GT(fp_t x, fp_t y) { return EQ(x,y) ? false : (x > y); }
	static bool GE(fp_t x, fp_t y) { return EQ(x,y) ? true  : (x > y); }

	static bool is_integer(fp_t x);	// is x fuzzily an integer?
	static int floor(fp_t x);	// round x fuzzily down to integer
	static int ceiling(fp_t x);	// round x fuzzily up to integer

private:
	// comparison tolerance
	// ... must be explicitly initialized when instantiating
	//     for a new <fp_t> type, see "fuzzy.cc" for details/examples
	static fp_t tolerance_;
	};

//******************************************************************************

//
// This template does machine-independent rounding of floating point
// values, parameterized by the floating point type.
//

// this template class has only static members
// ... it's a class, not a namespace, because we want to express the
//     semantics that the entire class is a single template, rather
//     than the individual members being conceptually-unrelated templates
template <typename fp_t>
class	round
	{
public:
	static int to_integer(fp_t x);	// round to nearest integer

	static int floor(fp_t x);	// round down to integer
	static int ceiling(fp_t x);	// round up to integer
	};

//******************************************************************************

	  }	// namespace jtutil

#endif	/* AHFINDERDIRECT__UTIL_HH */


