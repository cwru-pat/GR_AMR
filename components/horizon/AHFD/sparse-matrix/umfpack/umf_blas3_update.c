/* ========================================================================== */
/* === UMF_blas3_update ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

#include "umf_internal.h"

GLOBAL void UMF_blas3_update
(
    WorkType *Work
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Entry *A, *B, *C ;
    Int k, m, n, d ;

    /* ---------------------------------------------------------------------- */
    /* perform the matrix-matrix multiply */
    /* ---------------------------------------------------------------------- */

    DEBUG5 (("In UMF_blas3_update "ID" "ID" "ID"\n",
    	Work->fnpiv, Work->fnrows, Work->fncols)) ;

    k = Work->fnpiv ;
    Work->fnpiv = 0 ;	/* no more pivots in frontal working array */
    Work->fnzeros = 0 ;
    m = Work->fnrows ;
    n = Work->fncols ;
    if (k == 0 || m == 0 || n == 0)
    {
	return ;
    }

    d = Work->fnrows_max ;
    C = Work->Fx ;
    A =	C + (Work->fncols_max - k) * d ;
    B = C + d - k ;

    DEBUG5 (("DO RANK-NB UPDATE of frontal:\n")) ;
    DEBUG5 (("DGEMM : "ID" "ID" "ID"\n", k, m, n)) ;

    /* ---------------------------------------------------------------------- */
    /* C = C - A*B */
    /* ---------------------------------------------------------------------- */

    if (k == 1)
    {

#ifdef USE_NO_BLAS

	/* no BLAS available - use plain C code instead */
	Int i, j ;
	Entry b_sj, *c_ij, *a_is ;
	for (j = 0 ; j < n ; j++)
	{
	    b_sj = B [j*d] ;
	    if (IS_NONZERO (b_sj))
	    {
		c_ij = & C [j*d] ;
		a_is = & A [0] ;
		for (i = 0 ; i < m ; i++)
		{
		    /* (*c_ij++) -= b_sj * (*a_is++) ; */
		    MULT_SUB (*c_ij, b_sj, *a_is) ;
		    c_ij++ ;
		    a_is++ ;
		}
	    }
	}

#else

	BLAS_GER (m, n, A, B, C, d) ;

#endif /* USE_NO_BLAS */

    }
    else
    {

#ifdef USE_NO_BLAS

	/* no BLAS available - use plain C code instead */
	Int i, j, s ;
	Entry b_sj, *c_ij, *a_is ;
	for (j = 0 ; j < n ; j++)
	{
	    for (s = 0 ; s < k ; s++)
	    {
		b_sj = B [s+j*d] ;
		if (IS_NONZERO (b_sj))
		{
		    c_ij = & C [j*d] ;
		    a_is = & A [s*d] ;
		    for (i = 0 ; i < m ; i++)
		    {
			/* (*c_ij++) -= b_sj * (*a_is++) ; */
			MULT_SUB (*c_ij, b_sj, *a_is) ;
			c_ij++ ;
			a_is++ ;
		    }
		}
	    }
	}

#else

	BLAS_GEMM (m, n, k, A, B, C, d) ;

#endif /* USE_NO_BLAS */

    }

    DEBUG2 (("blas3 "ID" "ID" "ID"\n", k, Work->fnrows, Work->fncols)) ;
}

