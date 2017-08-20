// row_sparse_Jacobian.cc -- data structures for row-sparse Jacobian matrices
// $Header$
//
// <<<literal contents of "ilucg.h">>>
// data structures local to this file
// prototypes for functions local to this file
//
#ifdef HAVE_ROW_SPARSE_JACOBIAN
// row_sparse_Jacobian::row_sparse_Jacobian
// row_sparse_Jacobian::~row_sparse_Jacobian
// row_sparse_Jacobian::element
// row_sparse_Jacobian::zero_matrix
// row_sparse_Jacobian::set_element
// row_sparse_Jacobian::sum_into_element
/// row_sparse_Jacobian::find_element
/// row_sparse_Jacobian::insert_element
/// row_sparse_Jacobian::grow_arrays
/// row_sparse_Jacobian::sort_each_row_into_column_order
/// compare_matrix_elements
#ifdef DEBUG_ROW_SPARSE_JACOBIAN
/// row_sparse_Jacobian::check_and_print_data_structure
#endif
#endif
//
#ifdef HAVE_ROW_SPARSE_JACOBIAN__ILUCG
// row_sparse_Jacobian__ILUCG::row_sparse_Jacobian__ILUCG
// row_sparse_Jacobian__ILUCG::~row_sparse_Jacobian__ILUCG
// row_sparse_Jacobian__ILUCG::solve_linear_system
#endif
//
#ifdef HAVE_ROW_SPARSE_JACOBIAN__UMFPACK
// row_sparse_Jacobian__UMFPACK::row_sparse_Jacobian__UMFPACK
// row_sparse_Jacobian__UMFPACK::~row_sparse_Jacobian__UMFPACK
// row_sparse_Jacobian__UMFPACK::solve_linear_system
#endif
//

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "../jtutil/util_Table.h"
#include "../cctk.h"

using namespace SAMRAI;
#ifdef HAVE_ROW_SPARSE_JACOBIAN__ILUCG
//
// FIXME:
//   Cactus's CCTK_FCALL() isn't expanded in .h files (this is a bug),
//   so we include the contents of "../sparse-matrix/ilucg/ilucg.h"
//   inline below.  This is a doubleplusungood kludge! :( :(
//#include "../sparse-matrix/ilucg/ilucg.h"
//
#endif


#ifdef HAVE_ROW_SPARSE_JACOBIAN__UMFPACK
extern "C"
       {
  #include "../sparse-matrix/umfpack/umfpack.h"

  #ifdef FORTRAN_INTEGER_IS_INT
    #define umfpack_defaults		umfpack_di_defaults
    #define umfpack_symbolic		umfpack_di_symbolic
    #define umfpack_numeric		umfpack_di_numeric
    #define umfpack_wsolve		umfpack_di_wsolve
    #define umfpack_free_numeric	umfpack_di_free_numeric
    #define umfpack_free_symbolic	umfpack_di_free_symbolic
  #endif

  #ifdef FORTRAN_INTEGER_IS_LONG
    #define umfpack_defaults		umfpack_dl_defaults
    #define umfpack_symbolic		umfpack_dl_symbolic
    #define umfpack_numeric		umfpack_dl_numeric
    #define umfpack_wsolve		umfpack_dl_wsolve
    #define umfpack_free_numeric	umfpack_dl_free_numeric
    #define umfpack_free_symbolic	umfpack_dl_free_symbolic
  #endif
       }
#endif

#include "../jtutil/util.hh"
#include "../jtutil/array.hh"
#include "../jtutil/cpm_map.hh"
#include "../jtutil/linear_map.hh"

#include "../patch/coords.hh"
#include "../patch/grid.hh"
#include "../patch/fd_grid.hh"
#include "../patch/patch.hh"
#include "../patch/patch_edge.hh"
#include "../patch/patch_interp.hh"
#include "../patch/ghost_zone.hh"
#include "../patch/patch_system.hh"

#include "Jacobian.hh"
#include "row_sparse_Jacobian.hh"

// all the code in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// FIXME:
//   Cactus's CCTK_FCALL() isn't expanded in .h files (this is a bug),
//   so we include the contents of "../sparse-matrix/ilucg/ilucg.h"
//   inline here.  This is a doubleplusungood kludge! :( :(
//

//***** begin "ilucg.h" contents ******
/* ilucg.h -- C/C++ prototypes for [sd]ilucg() routines */
/* $Header$ */

/*
 * prerequisites:
 *	"cctk.h"
 *	"config.h"	// for "integer" = Fortran integer
 */

#ifdef __cplusplus
extern "C"
       {
#endif

/*
 * ***** ILUCG *****
 */
void CCTK_FCALL
  CCTK_FNAME(silucg)(const integer* N,
		     const integer IA[], const integer JA[], const float A[],
		     const float B[], float X[],
		     integer ITEMP[], float RTEMP[],
		     const float* eps, const integer* ITER,
		     integer* ISTATUS);
void CCTK_FCALL
  CCTK_FNAME(dilucg)(const integer* N,
		     const integer IA[], const integer JA[], const double A[],
		     const double B[], double X[],
		     integer ITEMP[], double RTEMP[],
		     const double* eps, const integer* ITER,
		     integer* ISTATUS);

#ifdef __cplusplus
       }	/* extern "C" */
#endif
//***** end "ilucg.h" contents ******

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// ***** data structures local to this file *****
//

// this represents a single element stored in the matrix for
// sort_row_into_column_order()  and  sort_row_into_column_order__cmp()
struct	matrix_element
	{
	int JA;
	fp  A;
	};

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// ***** prototypes for functions local to this file *****
//
namespace {
int compare_matrix_elements(const void* x, const void* y);
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN
//
// This function constructs a  row_sparse_Jacobian  object.
//
row_sparse_Jacobian::row_sparse_Jacobian(patch_system& ps, int IO_in,
					 bool print_msg_flag /* = false */)
	: Jacobian(ps), IO_(IO_in),
	#ifdef DEBUG_ROW_SPARSE_JACOBIAN
	  IOstr_((IO_in == 0) ? "c" : "f"),
	#endif
	  N_nonzeros_(0),
	  current_N_rows_(0),
	  N_nonzeros_allocated_(0),
	  IA_(new integer[N_rows_+1]),
	  JA_(NULL), A_(NULL)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "   row sparse matrix Jacobian [IO=%d] (%d rows)",
		   IO_, N_rows_);

zero_matrix();
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN */

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN
//
// This function destroys a  row_sparse_Jacobian  object.
//
row_sparse_Jacobian::~row_sparse_Jacobian()
{
delete[]  A_;
delete[] JA_;
delete[] IA_;
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN */

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN
//
// This function gets a matrix element.
//
fp row_sparse_Jacobian::element(int II, int JJ)
	const
{
const int posn = find_element(II,JJ);
return (posn >= 0) ? A_[posn] : 0.0;
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN */

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN
//
// This function zeros a  row_sparse_Jacobian  object.
//
void row_sparse_Jacobian::zero_matrix()
{
#ifdef ROW_SPARSE_JACOBIAN__LOG_MAIN_OPS
printf("row_sparse_Jacobian::zero_matrix()\n");
#endif

N_nonzeros_ = 0;
current_N_rows_= 0;
IA_[0] = IO_;
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN */

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN
//
// This function sets a matrix element to a specified value.
//
void row_sparse_Jacobian::set_element(int II, int JJ, fp value)
{
const int posn = find_element(II,JJ);
if (posn >= 0)
   then A_[posn] = value;
   else insert_element(II,JJ, value);
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN */

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN
//
// This function sums a value into a matrix element.
//
void row_sparse_Jacobian::sum_into_element(int II, int JJ, fp value)
{
const int posn = find_element(II,JJ);
if (posn >= 0)
   then A_[posn] += value;
   else insert_element(II,JJ, value);
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN */

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN
//
// This function searches our data structures for element (II,JJ).
// It returns the 0-origin position of the (II,JJ) element in A_ and JA_,
// or -1 if (II,JJ) is not found.
//
int row_sparse_Jacobian::find_element(int II, int JJ)
	const
{
if (II >= current_N_rows_)
   then return -1;			// this row not defined yet

const int start = IA_[II  ] - IO_;
const int stop  = IA_[II+1] - IO_;
	for (int posn = start ; posn < stop ; ++posn)
	{
	if (JA_[posn]-IO_ == JJ)
	   then return posn;				// found
	}

return -1;						// not found
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN */

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN
//
// This function inserts a new element in the matrix, either appending
// it to the current (last) row or starting a new row, and in either case
// growing the arrays if necessary.  It should only be called if II is
// the last row or one beyond, and JJ isn't already in row II.
//
int row_sparse_Jacobian::insert_element(int II, int JJ, fp value)
{
#ifdef ROW_SPARSE_JACOBIAN__LOG_MAIN_OPS
printf("row_sparse_Jacobian::insert_element(II=%d,JJ=%d, %s=%d,%d, value=%g)\n",
       II,JJ, IOstr_,II+IO_,JJ+IO_, value);
#endif


if (! ((II == current_N_rows_-1) || (II == current_N_rows_)) )
   then error_exit(PANIC_EXIT,
"***** row_sparse_Jacobian::insert_element(II=%d, JJ=%d, value=%g):\n"
"        attempt to insert element elsewhere than {last row, last row+1}!\n"
"        N_rows_=%d   current_N_rows_=%d   IO_=%d\n"
"        N_nonzeros_=%d   N_nonzeros_allocated_=%d\n"
		   ,
		   II, JJ, double(value),
		   N_rows_, current_N_rows_, IO_,
		   N_nonzeros_, N_nonzeros_allocated_);		/*NOTREACHED*/

// start a new row if necessary
if (II == current_N_rows_)
   then {
      #ifdef ROW_SPARSE_JACOBIAN__LOG_MAIN_OPS
	printf("   starting new row: current_N_rows_=%d", current_N_rows_);
      #endif
	assert(current_N_rows_ < N_rows_);
	IA_[current_N_rows_+1] = IA_[current_N_rows_];
	++current_N_rows_;
      #ifdef ROW_SPARSE_JACOBIAN__LOG_MAIN_OPS
	printf(" --> %d\n", current_N_rows_);
      #endif
      #ifdef ROW_SPARSE_JACOBIAN__CHECK_AT_ROW_START
	check_and_print_data_structure(false);
      #endif
      #ifdef ROW_SPARSE_JACOBIAN__PRINT_AT_ROW_START
	check_and_print_data_structure(true);
      #endif
	}

// insert into current row
assert(II == current_N_rows_-1);
if (IA_[II+1]-IO_ >= N_nonzeros_allocated_)
   then grow_arrays();
const int posn = IA_[II+1]-IO_;
assert(posn < N_nonzeros_allocated_);
JA_[posn] = JJ+IO_;
 A_[posn] = value;
++IA_[II+1];
++N_nonzeros_;
#ifdef ROW_SPARSE_JACOBIAN__LOG_MAIN_OPS
printf("   stored at posn=%d (new N_nonzeros_=%d)\n", posn, N_nonzeros_);
#endif

#ifdef ROW_SPARSE_JACOBIAN__CHECK_AFTER_INSERT
check_and_print_data_structure(false);
#endif
#ifdef ROW_SPARSE_JACOBIAN__PRINT_AFTER_INSERT
check_and_print_data_structure(true);
#endif

return posn;
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN */

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN
//
// This function grows the storage arrays to make room for more elements.
// The growth is geometric, so the amortized time to grow to size N is O(N).
//
void row_sparse_Jacobian::grow_arrays()
{
#ifdef ROW_SPARSE_JACOBIAN__LOG_MAIN_OPS
printf("   row_sparse_Jacobian::grow_arrays(): N_nonzeros_allocated=%d",
       N_nonzeros_allocated_);
#endif

N_nonzeros_allocated_ += base_growth_amount + (N_nonzeros_allocated_ >> 1);

#ifdef ROW_SPARSE_JACOBIAN__LOG_MAIN_OPS
printf(" --> %d\n", N_nonzeros_allocated_);
#endif

integer* const new_JA = new integer[N_nonzeros_allocated_];
fp     * const new_A  = new fp     [N_nonzeros_allocated_];
	for (int posn = 0 ; posn < N_nonzeros_ ; ++posn)
	{
	new_JA[posn] = JA_[posn];
	new_A [posn] = A_ [posn];
	}
delete[] A_ ;
delete[] JA_;
JA_ = new_JA;
A_  = new_A;
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN */

//******************************************************************************

//
// This function sorts each row's JA_[] and A_[] array elements
// so the columns (JA_[] values) are in increasing order.
// This doesn't change the abstract value of the matrix.
//
void row_sparse_Jacobian::sort_each_row_into_column_order()
{
//
// Implementation notes:
//
// The present implementation is quite inefficient: It copies the row's
// JA_[] and A_[] values into a contiguous buffer, calls qsort(3) on that
// buffer,
//	[the STL sort template would be a lot faster here,
//	but too many systems still have broken STLs -- sigh...]
// then copies the result back into JA_[] and A_[].
//

// buffer must be big enough to hold the largest row
int max_N_in_row = 0;
	  {
	for (int II = 0 ; II < N_rows_ ; ++II)
	{
        max_N_in_row = jtutil::max(max_N_in_row, IA_[II+1]-IA_[II]);
	}
	  }

// contiguous buffer for sorting
struct matrix_element *const buffer = new struct matrix_element[max_N_in_row];

	  {
	for (int II = 0 ; II < N_rows_ ; ++II)
	{
	const int N_in_row = IA_[II+1] - IA_[II];

	// copy this row's JA_[] and A_[] values to the buffer
	const int start = IA_[II]-IO_;
		for (int p = 0 ; p < N_in_row ; ++p)
		{
		const int posn = start + p;
		buffer[p].JA = JA_[posn];
		buffer[p].A  = A_ [posn];
		}

	// sort the buffer
	qsort(static_cast<void *>(buffer), N_in_row, sizeof(buffer[0]),
	      &compare_matrix_elements);

	// copy the buffer values back to this row's JA_[] and A_[]
		for (int p = 0 ; p < N_in_row ; ++p)
		{
		const int posn = start + p;
		JA_[posn] = buffer[p].JA;
		A_ [posn] = buffer[p].A;
		}
	}
	  }

delete[] buffer;
}

//**************************************

//
// This function is a qsort(3) comparison function for
// row_sparse_matrix::sort_each_row_into_column_order().
//
// Arguments:
// x,y = (void *) casts of pointers to  struct matrix_element  objects
//	 to be compared
//
// Results:
// This function returns < 0 if x < y, 0 if x == y, > 0 if x > y in the
// desired sort ordering.
namespace {
int compare_matrix_elements(const void* x, const void* y)
{
const struct matrix_element* const px
	= static_cast<const struct matrix_element*>(x);
const struct matrix_element* const py
	= static_cast<const struct matrix_element*>(y);

return px->JA - py->JA;
}
	  }

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN
#ifdef DEBUG_ROW_SPARSE_JACOBIAN
//
// This function sanity-checks and optionally prints the sparse-matrix
// data structures.
//
void row_sparse_Jacobian::check_and_print_data_structure(bool print_flag)
	const
{
if (print_flag)
   then {
	printf("--- begin Jacobian printout: "
	       "N_rows_=%d   current_N_rows_=%d   IO_=%d\n",
	       N_rows_, current_N_rows_, IO_);
	printf("                             "
	       "N_nonzeros_=%d   N_nonzeros_allocated_=%d\n",
		N_nonzeros_, N_nonzeros_allocated_);
	}

assert(N_rows_ >= 0);
assert(current_N_rows_ >= 0);
assert(current_N_rows_ <= N_rows_);
assert(N_nonzeros_ >= 0);
assert(N_nonzeros_allocated_ >= 0);
assert(N_nonzeros_ <= N_nonzeros_allocated_);

	for (int II = 0 ; II < current_N_rows_ ; ++II)
	{
	if (print_flag)
	   then printf("%sII=%d", IOstr_, II+IO_);

	const int start = IA_[II  ]-IO_;
	const int stop  = IA_[II+1]-IO_;
	if (print_flag)
	   then printf("\t%sposn=[%d,%d):", IOstr_, start+IO_, stop+IO_);
	if (II > 0)
	   then assert(IA_[II] >= IA_[II-1]);
	assert(start >= 0);
	assert(start <= N_nonzeros_allocated_);
	assert(stop >= 0);
	assert(stop <= N_nonzeros_allocated_);
	assert(stop >= start);

	if (print_flag)
	   then printf("\t%sJJ =", IOstr_);
		for (int posn = start ; posn < stop ; ++posn)
		{
		assert(posn >= 0);
		assert(posn < N_nonzeros_allocated_);
		const int gJJ = JA_[posn];
		if (print_flag)
		   then printf(" %d", gJJ);
		assert(gJJ >= IO_);
		assert(gJJ < N_rows_+IO_);
		}
	if (print_flag)
	   then printf("\n");
	}

if (print_flag)
   then printf("--- end Jacobian printout\n");
}
#endif	/* DEBUG_ROW_SPARSE_JACOBIAN */
#endif	/* HAVE_ROW_SPARSE_JACOBIAN */

//******************************************************************************
//******************************************************************************
//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN__ILUCG
//
// This function constructs a  row_sparse_Jacobian__ILUCG  object.
//
row_sparse_Jacobian__ILUCG::row_sparse_Jacobian__ILUCG
	(patch_system& ps,
	 bool print_msg_flag /* = false */)
	: row_sparse_Jacobian(ps, Fortran_index_origin,
			      print_msg_flag),
	  itemp_(NULL),
	  rtemp_(NULL)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      ILUCG linear-equations solver");
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN__ILUCG */

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN__ILUCG
//
// This function destroys a  row_sparse_Jacobian__ILUCG  object.
//
row_sparse_Jacobian__ILUCG::~row_sparse_Jacobian__ILUCG()
{
delete[] rtemp_;
delete[] itemp_;
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN__ILUCG */

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN__ILUCG
//
// This function solves the linear system J.x = rhs, with rhs and x
// being nominal-grid gridfns, using an ILUCG subroutine originally
// provided to me in 1985 (!) by Tom Nicol of the UBC Computing Center.
//
// It returns -1.0 (no condition number estimate).
//
// FIXME:
// It would be more efficient to adjust the ILUCG error tolerance
// based on the current Newton-iteration error level (in ../driver/Newton.cc).
// That is, at early stages of the Newton iteration there's no need to
// solve this linear system to high accuracy.
//
fp row_sparse_Jacobian__ILUCG::solve_linear_system
	(int rhs_gfn, int x_gfn,
	 const struct linear_solver_pars& pars,
	 bool print_msg_flag)
{
assert(IO_ == Fortran_index_origin);	// ILUCG uses Fortran indices
assert(current_N_rows_ == N_rows_);	// matrix must be fully defined

//
// if this is our first call, allocate the scratch arrays
//
if (itemp_ == NULL)
   then {
	if (print_msg_flag)
	   then {
		CCTK_VInfo(CCTK_THORNSTRING,
			   "row_sparse_Jacobian__ILUCG::solve_linear_system()");
		CCTK_VInfo(CCTK_THORNSTRING,
		   "   N_rows_=%d N_nonzeros_=%d N_nonzeros_allocated_=%d",
			   N_rows_, N_nonzeros_, N_nonzeros_allocated_);
		}
	itemp_ = new integer[3*N_rows_ + 3*N_nonzeros_ + 2];
	rtemp_ = new fp     [4*N_rows_ +   N_nonzeros_    ];
	}

//
// set up the ILUCG subroutine
//

// initial guess = all zeros
fp *x = ps_.gridfn_data(x_gfn);
	for (int II = 0 ; II < N_rows_ ; ++II)
	{
	x[II] = 0.0;
	}

const integer N = N_rows_;
const fp *rhs = ps_.gridfn_data(rhs_gfn);
const fp      eps            = pars.ILUCG_pars.error_tolerance;
const integer max_iterations = pars.ILUCG_pars.limit_CG_iterations
			       ? N_rows_ : 0;
integer istatus;

// the actual linear solution
#if   defined(FP_IS_FLOAT)
  CCTK_FNAME(silucg)(&N,
		     IA_, JA_, A_,
		     rhs, x,
		     itemp_, rtemp_,
		     &eps, &max_iterations,
		     &istatus);
#elif defined(FP_IS_DOUBLE)
  CCTK_FNAME(dilucg)(&N,
		     IA_, JA_, A_,
		     rhs, x,
		     itemp_, rtemp_,
		     &eps, &max_iterations,
		     &istatus);
#else
  #error "don't know fp datatype!"
#endif
if (istatus < 0)
   then error_exit(ERROR_EXIT,
"***** row_sparse_Jacobian__ILUCG::solve_linear_system(rhs_gfn=%d, x_gfn=%d):\n"
"        error return from [sd]ilucg() routine!\n"
"        istatus=%d < 0 ==> bad matrix structure, eg. zero diagonal element!\n"
		   ,
		   rhs_gfn, x_gfn,
		   int(istatus));				/*NOTREACHED*/

//
// If istatus == 0,
//    we didn't reach the convergence level within the specified
//    maximum number of conjugate gradient (CG) iterations. :(
//
//    In practice, however, most of the time we don't need an
//    uber-accurate solution, so it's reasonable to just continue here
//    with whatever approximate solution the CG iteration reached.
//
//    If the Newton iteration converges, we're fine (if |Theta| is
//    small then we have a good horizon position regardless of how we
//    got there), so the only danger in continuing with an approximate
//    linear-system solution here is that we may slow or even prevent
//    the convergence of the Newton iteration.  In practice, this
//    doesn't seem to be a serious problem...
//
const int actual_N_iterations = (istatus == 0) ? max_iterations : istatus;
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "   %d CG iteration%s%s",
		   actual_N_iterations,
		   ((actual_N_iterations == 1) ? "" : "s"),
		   ((istatus == 0) ? " (no convergence ==> continuing)"
				   : " (converged ok)"));

return -1.0;			// no condition number estimate available
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN__ILUCG */

//******************************************************************************
//******************************************************************************
//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN__UMFPACK
//
// This function constructs a  row_sparse_Jacobian__UMFPACK  object.
//
row_sparse_Jacobian__UMFPACK::row_sparse_Jacobian__UMFPACK
	(patch_system& ps,
	 bool print_msg_flag /* = false */)
	: row_sparse_Jacobian(ps, C_index_origin,
			      print_msg_flag),
	  Control_(new double[UMFPACK_CONTROL]),
	  Info_   (new double[UMFPACK_INFO   ]),
	  solve_workspace_integer_(NULL),
	  solve_workspace_double_ (NULL),
	  Symbolic_(NULL),
	  Numeric_(NULL)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      UMFPACK linear-equations solver");

umfpack_defaults(Control_);
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN__UMFPACK */

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN__UMFPACK
//
// This function destroys a  row_sparse_Jacobian__UMFPACK  object.
//
row_sparse_Jacobian__UMFPACK::~row_sparse_Jacobian__UMFPACK()
{
if (Numeric_ != NULL)
   then umfpack_free_numeric (&Numeric_ );
if (Symbolic_ != NULL)
   then umfpack_free_symbolic(&Symbolic_);
delete[] solve_workspace_double_;
delete[] solve_workspace_integer_;
delete[] Info_;
delete[] Control_;
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN__UMFPACK */

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN__UMFPACK
//
// This function solves the linear system J.x = rhs, with rhs and x
// being nominal-grid gridfns, using the UMFPACK sparse-LU-decomposition
// package.
//
// Because the UMFPACK condition number estimator doesn't seem to be
// reliable, this function conditionally prints it, but returns -1.0
// to indicate that no (reliable) condition number estimate is available.
//
fp row_sparse_Jacobian__UMFPACK::solve_linear_system
	(int rhs_gfn, int x_gfn,
	 const struct linear_solver_pars& pars,
	 bool print_msg_flag)
{
int status;

assert(IO_ == C_index_origin);		// UMFPACK uses C indices
assert(current_N_rows_ == N_rows_);	// matrix must be fully defined

if (print_msg_flag)
   then {
	CCTK_VInfo(CCTK_THORNSTRING,
		   "row_sparse_Jacobian__UMFPACK::solve_linear_system()");
	CCTK_VInfo(CCTK_THORNSTRING,
		   "   N_rows_=%d N_nonzeros_=%d N_nonzeros_allocated_=%d",
		   N_rows_, N_nonzeros_, N_nonzeros_allocated_);
	}

//
// we set up the matrix in row-sparse format, with the entries in each
// column unsorted, but UMFPACK wants column-sparse format, with the
// entries in each column sorted in increasing order of their rows
// ==> sort the entries in each row,
//     then use UMFPACK to solve the transposed system
//
sort_each_row_into_column_order();

//
// if we haven't done so already, symbolically factor the matrix
//
if (Symbolic_ == NULL)
   then {
	if (print_msg_flag)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "   UMFPACK symbolic factorization");
	status = umfpack_symbolic(N_rows_, N_rows_,	// matrix size
				  IA_, JA_,		// sparsity structure
				  &Symbolic_,	// umfpack sets Symbolic to
						// point to the new object
				  Control_, Info_);
	if (status != UMFPACK_OK)
	   then error_exit(ERROR_EXIT,
"***** row_sparse_Jacobian__UMFPACK::solve_linear_system():\n"
"        error return status=%d from umfpack_symbolic() routine\n"
		   ,
			   status);				/*NOTREACHED*/
	}

//
// do the sparse LU decomposition
//
assert(Symbolic_ != NULL);
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "   UMFPACK sparse LU decomposition");
status = umfpack_numeric(IA_, JA_, A_,		// input matrix
			 Symbolic_,
			 &Numeric_,
			 Control_, Info_);
if (status != UMFPACK_OK)
   then error_exit(ERROR_EXIT,
"***** row_sparse_Jacobian__UMFPACK::solve_linear_system():\n"
"        error return status=%d from umfpack_numeric() routine\n"
		   ,
		   status);					/*NOTREACHED*/
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      UMFPACK rcond=%.1e (unreliable ==> ignoring)",
		   double(Info_[UMFPACK_RCOND]));
const fp rcond = -1.0;


//
// set parameters for linear solve
//
if (pars.UMFPACK_pars.N_II_iterations >= 0)
   then Control_[UMFPACK_IRSTEP] = pars.UMFPACK_pars.N_II_iterations;
   else {
	// accept the UMFPACK default (which we've already set)
	// ==> no-op here
	}
const bool II_flag = (Control_[UMFPACK_IRSTEP] > 0);


//
// if we haven't done so already, allocate workspace
// (how much depends on whether or not UMFPACK does iterative improvement)
//
if (solve_workspace_integer_ == NULL)
   then solve_workspace_integer_ = new integer[N_rows_];
if (solve_workspace_double_ == NULL)
   then solve_workspace_double_  = new double [II_flag ? 5*N_rows_ : N_rows_];


//
// solve the transposed linear system
// 
const double* rhs = ps_.gridfn_data(rhs_gfn);
      double*   x = ps_.gridfn_data(  x_gfn);
if (print_msg_flag)
   then {
	CCTK_VInfo(CCTK_THORNSTRING,
		   "   UMFPACK solution using LU decomposition");
	if (II_flag)
	   then CCTK_VInfo(CCTK_THORNSTRING,
		   "                          and iterative improvement");
	}
status = umfpack_wsolve(UMFPACK_At,		// solve transposed system
			IA_, JA_, A_,		// input matrix
			x,			// solution vector
			rhs,			// rhs vector
			Numeric_,		// LU decomposition object
			Control_, Info_,
			solve_workspace_integer_,
			solve_workspace_double_);
if (status != UMFPACK_OK)
   then error_exit(ERROR_EXIT,
"***** row_sparse_Jacobian__UMFPACK::solve_linear_system():\n"
"        error return status=%d from umfpack_wsolve() routine\n"
		   ,
		   status);					/*NOTREACHED*/

if (Numeric_ != NULL)
   then umfpack_free_numeric (&Numeric_ );

return rcond;
}
#endif	/* HAVE_ROW_SPARSE_JACOBIAN__UMFPACK */

//******************************************************************************
//******************************************************************************
//******************************************************************************

	  }	// namespace AHFinderDirect
