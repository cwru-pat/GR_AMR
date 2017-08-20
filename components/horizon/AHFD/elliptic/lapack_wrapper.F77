c lapack_wrapper.f -- wrapper routines for LAPACK [sd]gecon()
c $Header$

c
c These subroutines are wrappers around the LAPACK [sd]gecon() subroutines.
c These subroutines take only integer/real/double precision arguments,
c avoiding problems with C/C++ --> Fortran passing of the character string
c arguments used by [sd]gecon().
c
c Arguments:
c norm_int = (in) 0 ==> infinity-norm
c                 1 ==> 1-norm
c

c
c Note that the Compaq f90 compiler complains about empty (or all-comment)
c files, so we still have to define an empty subroutine even if it is never
c called.  :( Oh well, memory is cheap...
c

#include "../AHFD_macros.h"

	subroutine sgecon_wrapper(norm_int,
     $	                          N, A, LDA, anorm, rcond,
     $	                          WORK, IWORK, info)
	integer norm_int
	integer N, LDA
	real A(LDA,N)
	real anorm, rcond
	real WORK(*)
	integer iwork(*)
	integer info
#ifdef HAVE_DENSE_JACOBIAN__LAPACK
#ifdef FP_IS_FLOAT
	if	(norm_int .eq. 0) then
		call sgecon('I', N,A,LDA, anorm,rcond, WORK,IWORK, info)
	else if (norm_int .eq. 1) then
		call sgecon('1', N,A,LDA, anorm,rcond, WORK,IWORK, info)
	else
		info = -1;
	end if
#endif /* FP_IS_FLOAT */
#endif	/* HAVE_DENSE_JACOBIAN__LAPACK */
	return
	end

	subroutine dgecon_wrapper(norm_int,
     $	                          N, A, LDA, anorm, rcond,
     $	                          WORK, IWORK, info)
	integer norm_int
	integer N, LDA
	double precision A(LDA,N)
	double precision anorm, rcond
	double precision WORK(*)
	integer iwork(*)
	integer info
#ifdef HAVE_DENSE_JACOBIAN__LAPACK
#ifdef FP_IS_DOUBLE
	if	(norm_int .eq. 0) then
		call dgecon('I', N,A,LDA, anorm,rcond, WORK,IWORK, info)
	else if (norm_int .eq. 1) then
		call dgecon('1', N,A,LDA, anorm,rcond, WORK,IWORK, info)
	else
		info = -1;
	end if
#endif /* FP_IS_DOUBLE */
#endif	/* HAVE_DENSE_JACOBIAN__LAPACK */
	return
	end
