// dense_Jacobian.cc -- dense-matrix Jacobian
// $Header$

//
// <<<literal contents of "lapack.hh">>>
//
#ifdef HAVE_DENSE_JACOBIAN
// dense_Jacobian::dense_Jacobian
// dense_Jacobian::zero_matrix
#endif
//
#ifdef HAVE_DENSE_JACOBIAN__LAPACK
// dense_Jacobian__LAPACK::dense_Jacobian__LAPACK
// dense_Jacobian__LAPACK::~dense_Jacobian__LAPACK
// dense_Jacobian__LAPACK::solve_linear_system
#endif
//

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "../jtutil/util_Table.h"
#include "../AHFD_macros.h"

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
#include "dense_Jacobian.hh"

// all the code in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// FIXME:
//   Cactus's CCTK_FCALL() isn't expanded in .h files (this is a bug),
//   so we include the contents of "lapack.h" inline here.  This is a
//   kludge! :( :(
//#include "lapack.h"
//

//***** begin "lapack.h" contents ******
/* lapack.h -- C/C++ prototypes for (some) BLAS+LAPACK+wrapper routines */
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
 * ***** BLAS *****
 */
integer CCTK_FCALL
  CCTK_FNAME(isamax)(const integer* N, const float SX[], const integer* incx);
integer CCTK_FCALL
  CCTK_FNAME(idamax)(const integer* N, const double DX[], const integer* incx);

/*
 * ***** LAPACK *****
 */
void CCTK_FCALL
  CCTK_FNAME(sgesv)(const integer* N, const integer* NRHS,
		    float A[], const int* LDA,
		    integer IPIV[],
		    float B[], const integer* LDB, integer* info);
void CCTK_FCALL
  CCTK_FNAME(dgesv)(const integer* N, const integer* NRHS,
		    double A[], const int* LDA,
		    integer IPIV[],
		    double B[], const integer* LDB, integer* info);

/*
 * ***** wrappers (for passing character-string args) *****
 */
/* norm_int = 0 for infinity-norm, 1 for 1-norm */
void CCTK_FCALL
  CCTK_FNAME(sgecon_wrapper)(const integer* norm_int,
			     const integer* N,
			     float A[], const integer* LDA,
			     const float* anorm, float* rcond,
			     float WORK[], integer IWORK[],
			     integer* info);
void CCTK_FCALL
  CCTK_FNAME(dgecon_wrapper)(const integer* norm_int,
			     const integer* N,
			     double A[], const integer* LDA,
			     const double* anorm, double* rcond,
			     double WORK[], integer IWORK[],
			     integer* info);

#ifdef __cplusplus
       }	/* extern "C" */
#endif
//***** end "lapack.h" contents ********

//******************************************************************************
//******************************************************************************
//******************************************************************************

#ifdef HAVE_DENSE_JACOBIAN
//
// This function constructs a  dense_Jacobian  object.
//
dense_Jacobian::dense_Jacobian(patch_system& ps,
			       bool print_msg_flag /* = false */)
	: Jacobian(ps),
	  matrix_(0,N_rows_-1, 0,N_rows_-1)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "   dense Jacobian matrix (%d rows)",
		   N_rows_);
}
#endif	/* HAVE_DENSE_JACOBIAN */

//******************************************************************************

#ifdef HAVE_DENSE_JACOBIAN
//
// This function zeros a  dense_Jacobian  object.
// Bugs: it could be made more efficient...
//
void dense_Jacobian::zero_matrix()
{
	for (int JJ = 0 ; JJ < N_rows_ ; ++JJ)
	{
		for (int II = 0 ; II < N_rows_ ; ++II)
		{
		set_element(II,JJ, 0.0);
		}
	}
}
#endif	/* HAVE_DENSE_JACOBIAN */

//******************************************************************************
//******************************************************************************
//******************************************************************************

#ifdef HAVE_DENSE_JACOBIAN__LAPACK
//
// This function constructs a  dense_Jacobian__LAPACK  object.
//
dense_Jacobian__LAPACK::dense_Jacobian__LAPACK
	(patch_system& ps,
	 bool print_msg_flag /* = false */)
	: dense_Jacobian(ps, print_msg_flag),
	  pivot_(new integer[  N_rows_]),
	  iwork_(new integer[  N_rows_]),
	  rwork_(new fp     [4*N_rows_])
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      LAPACK linear-equations solver");
}
#endif	/* HAVE_DENSE_JACOBIAN__LAPACK */

//******************************************************************************

#ifdef HAVE_DENSE_JACOBIAN__LAPACK
//
// This function destroys a  dense_Jacobian__LAPACK  object.
//
dense_Jacobian__LAPACK::~dense_Jacobian__LAPACK()
{
delete[] rwork_;
delete[] iwork_;
delete[] pivot_;
}
#endif	/* HAVE_DENSE_JACOBIAN__LAPACK */

//******************************************************************************

#ifdef HAVE_DENSE_JACOBIAN__LAPACK
//
// This function solves the linear system J.x = rhs, with rhs and x
// being nominal-grid gridfns, using LAPACK LU-decomposition routines.
//
// It returns the (infinity-norm) estimated reciprocal condition number
// of the linear system.
//
fp dense_Jacobian__LAPACK::solve_linear_system
	(int rhs_gfn, int x_gfn,
	 const struct linear_solver_pars& pars,
	 bool print_msg_flag)
{
const fp *rhs = ps_.gridfn_data(rhs_gfn);
fp       *x   = ps_.gridfn_data(x_gfn);
fp       *A   = matrix_.data_array();

const integer infinity_norm_flag = 0;
const integer stride = 1;
const integer NRHS = 1;
const integer N = N_rows_;
const integer N2 = N_rows_ * N_rows_;

// compute the infinity-norm of the matrix A
// ... max_posn = 1-origin index of A[] element with largest absolute value
#if   defined(FP_IS_FLOAT)
  const integer max_posn = CCTK_FNAME(isamax)(&N2, A, &stride);
#elif defined(FP_IS_DOUBLE)
  const integer max_posn = CCTK_FNAME(idamax)(&N2, A, &stride);
#else
  #error "don't know fp datatype!"
#endif
const fp A_infnorm = jtutil::abs(A[max_posn-1]);

// LU decompose and solve the linear system
//
// ... [sd]gesv() use an "in out" design, where the same argument
//     is used for both rhs and x ==> we must first copy rhs to x
//
	for (int II = 0 ; II < N_rows_ ; ++II)
	{
	x[II] = rhs[II];
	}
integer info;
#if   defined(FP_IS_FLOAT)
  CCTK_FNAME(sgesv)(&N, &NRHS, A, &N, pivot_, x, &N, &info);
#elif defined(FP_IS_DOUBLE)
  CCTK_FNAME(dgesv)(&N, &NRHS, A, &N, pivot_, x, &N, &info);
#else
  #error "don't know fp datatype!"
#endif

if (info < 0)
   then error_exit(ERROR_EXIT,
"\n"
"***** dense_Jacobian__LAPACK::solve_linear_system(rhs_gfn=%d, x_gfn=%d):\n"
"        error return (bad argument) info=%d from [sd]gesv() LAPACK routine!"
		   ,
		   rhs_gfn, x_gfn,
		   int(info));					/*NOTREACHED*/

if (info > 0)
   then return 0.0;					// *** ERROR RETURN ***
							// *** (singular matrix)

// estimate infinity-norm condition number
fp rcond;
#if   defined(FP_IS_FLOAT)
  CCTK_FNAME(sgecon_wrapper)(&infinity_norm_flag,
			     &N, A, &N, &A_infnorm, &rcond,
			     rwork_, iwork_, &info);
#elif defined(FP_IS_DOUBLE)
  CCTK_FNAME(dgecon_wrapper)(&infinity_norm_flag,
			     &N, A, &N, &A_infnorm, &rcond,
			     rwork_, iwork_, &info);
#else
  #error "don't know fp datatype!"
#endif

return rcond;
}
#endif	/* HAVE_DENSE_JACOBIAN__LAPACK */

//******************************************************************************
//******************************************************************************
//******************************************************************************

	  }	// namespace AHFinderDirect
