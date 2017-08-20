// array.cc -- array template classes
// $Header$
//
// jtutil::array1d::array1d - 1D array template constructor
// jtutil::array1d::~array1d - 1D array template destructor
//
// jtutil::array2d::array2d - 2D array template constructor
// jtutil::array2d::~array2d - 2D array template destructor
//
// jtutil::array3d::array3d - 3D array template constructor
// jtutil::array3d::~array3d - 3D array template destructor
//
#ifdef NOT_USED
// jtutil::array4d::array4d - 4D array template constructor
// jtutil::array4d::~array4d - 4D array template destructor
#endif
//
// ***** template instantiations *****
//

#include <stddef.h>	// NULL
#include <stdlib.h>	// size_t
#include <cmath>

// we want to instantiate templates with CCTK_* types
#include "../cctk.h"


#include "util.hh"
#include "array.hh"

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function constructs an  array1d  object.
//
namespace jtutil
	  {
template <typename T>
array1d<T>::array1d(int min_i_in, int max_i_in,
		    T *array_in /* = NULL */,
		    int stride_i_in /* = 0 */)
	: array_(array_in),
	  offset_(0),		// temp value, changed below
	  stride_i_(stride_i_in),
	  min_i_(min_i_in), max_i_(max_i_in),
	  we_own_array_(array_in == NULL)
{
if (stride_i_ == 0)
   then stride_i_ = 1;

// must use unchecked subscripting here since setup isn't done yet
offset_ = - subscript_unchecked(min_i_);		// RHS uses offset_ = 0
TBOX_ASSERT( subscript_unchecked(min_i_) == 0 );
max_subscript_ = subscript_unchecked(max_i_);

if (we_own_array_)
   then {
	// allocate it
	const int N_allocate = N_i();
	array_ = new T[N_allocate];
	}

// explicitly initialize array (new[] *doesn't* do this automagically)
	for (int i = min_i() ; i <= max_i() ; ++i)
	{
	operator()(i) = T(0);
	}
}
	  }	// namespace jtutil::

//******************************************************************************

//
// This function destroys an  array1d  object.
//
namespace jtutil
	  {
template <typename T>
array1d<T>::~array1d()
{
if (we_own_array_)
   then delete[] array_;
}
	  }	// namespace jtutil::

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function constructs an  array2d  object.
//
namespace jtutil
	  {
template <typename T>
array2d<T>::array2d(int min_i_in, int max_i_in,
		    int min_j_in, int max_j_in,
		    T *array_in /* = NULL */,
		    int stride_i_in /* = 0 */, int stride_j_in /* = 0 */)
	: array_(array_in),
	  offset_(0),		// temp value, changed below
	  stride_i_(stride_i_in), stride_j_(stride_j_in),
	  min_i_(min_i_in), max_i_(max_i_in),
	  min_j_(min_j_in), max_j_(max_j_in),
	  we_own_array_(array_in == NULL)
{
if (stride_j_ == 0)
   then stride_j_ = 1;
if (stride_i_ == 0)
   then stride_i_ = N_j();

// must use unchecked subscripting here since setup isn't done yet
offset_ = - subscript_unchecked(min_i_,min_j_);		// RHS uses offset_ = 0
assert( subscript_unchecked(min_i_,min_j_) == 0 );
max_subscript_ = subscript_unchecked(max_i_,max_j_);

if (we_own_array_)
   then {
	// allocate it
	const int N_allocate = N_i() * N_j();
	array_ = new T[N_allocate];
	}

// explicitly initialize array (new[] *doesn't* do this automagically)
	for (int i = min_i() ; i <= max_i() ; ++i)
	{
	for (int j = min_j() ; j <= max_j() ; ++j)
	{
	operator()(i,j) = T(0);
	}
	}
}
	  }	// namespace jtutil::

//******************************************************************************

//
// This function destroys an  array2d  object.
//
namespace jtutil
	  {
template <typename T>
array2d<T>::~array2d()
{
if (we_own_array_)
   then delete[] array_;
}
	  }	// namespace jtutil::

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function constructs an  array3d  object.
//
namespace jtutil
	  {
template <typename T>
array3d<T>::array3d(int min_i_in, int max_i_in,
		    int min_j_in, int max_j_in,
		    int min_k_in, int max_k_in,
		    T *array_in /* = NULL */,
		    int stride_i_in /* = 0 */, int stride_j_in /* = 0 */,
		    int stride_k_in /* = 0 */)
	: array_(array_in),
	  offset_(0),		// temp value, changed below
	  stride_i_(stride_i_in), stride_j_(stride_j_in),
	  stride_k_(stride_k_in),
	  min_i_(min_i_in), max_i_(max_i_in),
	  min_j_(min_j_in), max_j_(max_j_in),
	  min_k_(min_k_in), max_k_(max_k_in),
	  we_own_array_(array_in == NULL)
{
if (stride_k_ == 0)
   then stride_k_ = 1;
if (stride_j_ == 0)
   then stride_j_ = N_k();
if (stride_i_ == 0)
   then stride_i_ = N_j()*N_k();

// must use unchecked subscripting here since setup isn't done yet
offset_ = - subscript_unchecked(min_i_,min_j_,min_k_);	// RHS uses offset_ = 0
assert( subscript_unchecked(min_i_,min_j_,min_k_) == 0 );
max_subscript_ = subscript_unchecked(max_i_,max_j_,max_k_);

if (we_own_array_)
   then {
	// allocate it
	const int N_allocate = N_i() * N_j() * N_k();
	array_ = new T[N_allocate];
	}

// explicitly initialize array (new[] *doesn't* do this automagically)
	for (int i = min_i() ; i <= max_i() ; ++i)
	{
	for (int j = min_j() ; j <= max_j() ; ++j)
	{
	for (int k = min_k() ; k <= max_k() ; ++k)
	{
	operator()(i,j,k) = T(0);
	}
	}
	}
}
	  }	// namespace jtutil::

//******************************************************************************

//
// This function destroys an  array3d  object.
//
namespace jtutil
	  {
template <typename T>
array3d<T>::~array3d()
{
if (we_own_array_)
   then delete[] array_;
}
	  }	// namespace jtutil::

//******************************************************************************
//******************************************************************************
//******************************************************************************

#ifdef NOT_USED
//
// This function constructs an  array4d  object.
//
namespace jtutil
	  {
template <typename T>
array4d<T>::array4d(int min_i_in, int max_i_in,
		    int min_j_in, int max_j_in,
		    int min_k_in, int max_k_in,
		    int min_l_in, int max_l_in,
		    T *array_in /* = NULL */,
		    int stride_i_in /* = 0 */, int stride_j_in /* = 0 */,
		    int stride_k_in /* = 0 */, int stride_l_in /* = 0 */)
	: array_(array_in),
	  offset_(0),		// temp value, changed below
	  stride_i_(stride_i_in), stride_j_(stride_j_in),
	  stride_k_(stride_k_in), stride_l_(stride_l_in),
	  min_i_(min_i_in), max_i_(max_i_in),
	  min_j_(min_j_in), max_j_(max_j_in),
	  min_k_(min_k_in), max_k_(max_k_in),
	  min_l_(min_l_in), max_l_(max_l_in),
	  we_own_array_(array_in == NULL)
{
if (stride_l_ == 0)
   then stride_l_ = 1;
if (stride_k_ == 0)
   then stride_k_ = N_l();
if (stride_j_ == 0)
   then stride_j_ = N_k()*N_l();
if (stride_i_ == 0)
   then stride_i_ = N_j()*N_k()*N_l();

// must use unchecked subscripting here since setup isn't done yet
offset_ = - subscript_unchecked(min_i_,min_j_,
				min_k_,min_l_);		// RHS uses offset_ = 0
assert( subscript_unchecked(min_i_,min_j_,min_k_,min_l_) == 0 );
max_subscript_ = subscript_unchecked(max_i_,max_j_,max_k_,max_l_);

if (we_own_array_)
   then {
	// allocate it
	const int N_allocate = N_i() * N_j() * N_k() * N_l();
	array_ = new T[N_allocate];
	}

// explicitly initialize array (new[] *doesn't* do this automagically)
	for (int i = min_i() ; i <= max_i() ; ++i)
	{
	for (int j = min_j() ; j <= max_j() ; ++j)
	{
	for (int k = min_k() ; k <= max_k() ; ++k)
	{
	for (int l = min_l() ; l <= max_l() ; ++l)
	{
	operator()(i,j,k,l) = T(0);
	}
	}
	}
	}
}
	  }	// namespace jtutil::
#endif	/* NOT_USED */

//******************************************************************************

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// ***** template instantiations *****
//

template class jtutil::array1d<int>;
template class jtutil::array1d<long>;

// FIXME: we shouldn't have to instantiate these both, the const one
//	  is actually trivially derivable from the non-const one. :(
template class jtutil::array1d<      void *>;
template class jtutil::array1d<const void *>;

// full-fledged Cactus thorn
template class jtutil::array1d<CCTK_REAL>;

template class jtutil::array2d<CCTK_INT>;
template class jtutil::array2d<CCTK_REAL>;

template class jtutil::array3d<CCTK_REAL>;
