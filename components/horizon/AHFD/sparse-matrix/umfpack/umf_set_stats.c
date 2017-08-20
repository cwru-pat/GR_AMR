/* ========================================================================== */
/* === UMF_set_stats ======================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Sets statistics in Info array.  Calculates everything in double precision,
    rather than Int or size_t, so that usage estimates can be computed even if
    the problem is so large that it would cause integer overflow.

    This routine has many double relop's, but the NaN case is ignored.
*/

#include "umf_internal.h"
#include "umf_symbolic_usage.h"

GLOBAL void UMF_set_stats
(
    double Info [ ],
    SymbolicType *Symbolic,
    double max_usage,		/* peak size of Numeric->Memory, in Units */
    double num_mem_size,	/* final size of Numeric->Memory, in Units */
    double flops,		/* "true flops" */
    double lnz,			/* nz in L */
    double unz,			/* nz in U */
    double maxfrsize,		/* largest front size */
    double ulen,		/* size of Numeric->Upattern */
    double npiv,		/* number of pivots found */
    Int what			/* ESTIMATE or ACTUAL */
)
{

    double sym_size, work_usage, nn, n_row, n_col, n_inner, num_On_size1,
	num_On_size2, num_usage, maxncols, maxnrows ;

    n_col = Symbolic->n_col ;
    n_row = Symbolic->n_row ;
    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;
    maxncols = Symbolic->maxncols ;
    maxnrows = Symbolic->maxnrows ;

    /* final Symbolic object size */
    sym_size = UMF_symbolic_usage (Symbolic->n_row, Symbolic->n_col,
	Symbolic->nchains, Symbolic->nfr) ;

    /* size of O(n) part of Numeric object during factorization, */
    /* except Numeric->Memory and Numeric->Upattern */
    num_On_size1 =
	DUNITS (NumericType, 1)		/* Numeric structure */
	+ DUNITS (Entry, n_inner+1)	/* D */
	+ 4 * DUNITS (Int, n_row+1)	/* Rperm, Lpos, Uilen, Uip */
	+ 4 * DUNITS (Int, n_col+1) ;	/* Cperm, Upos, Lilen, Lip */

    /* size of O(n) part of Numeric object after factorization, */
    /* except Numeric->Memory and Numeric->Upattern */
    num_On_size2 =
	DUNITS (NumericType, 1)		/* Numeric structure */
	+ DUNITS (Entry, n_inner+1)	/* D */
	+ DUNITS (Int, n_row+1)		/* Rperm */
	+ DUNITS (Int, n_col+1)		/* Cperm */
	+ 4 * DUNITS (Int, npiv+1) ;	/* Lpos, Uilen, Uip, Upos, Lilen, Lip */

    /* peak size of Numeric->Memory */
    Info [UMFPACK_VARIABLE_PEAK + what] = max_usage ;

    /* final size of Numeric->Memory */
    Info [UMFPACK_VARIABLE_FINAL + what] = num_mem_size ;

    /* final size of Numeric object, including Numeric->Memory and ->Upattern */
    Info [UMFPACK_NUMERIC_SIZE + what] =
	num_On_size2
	+ num_mem_size		/* final Numeric->Memory size */
	+ DUNITS (Int, ulen) ;	/* Numeric->Upattern (from Work->Upattern) */

    /* largest front size (Work->Fx size, or actual size used) */
    Info [UMFPACK_MAX_FRONT_SIZE + what] = maxfrsize ;

    /* UMF_kernel work usage */
    work_usage =
	/* Work-> arrays */
	DUNITS (Entry, Symbolic->maxfrsize)	/* Fx */
	+ 2 * DUNITS (Int, n_row+1)		/* Frpos, Lpattern */
	+ 2 * DUNITS (Int, n_col+1)		/* Fcpos, Upattern */
	+ DUNITS (Int, n_col+n_inner+1)		/* E */
	+ DUNITS (Int, nn + 1)			/* Wp */
	+ 3 * DUNITS (Int, maxncols + 1)	/* Fcols, Wio, Woi */
	+ 2 * DUNITS (Int, maxnrows + 1)	/* Frows, Wm */
	+ DUNITS (Int, MAX (maxnrows, maxncols) + 1)	/* Woo */
	+ DUNITS (Int, Symbolic->nfr + 1) ;	/* Front_new1strow */

    /* Peak memory for just UMFPACK_numeric.  This excludes Numeric->Upattern */
    /* since it includes the equivalenced Work->Upattern array. */
    num_usage = sym_size + num_On_size1 + work_usage + max_usage ;

    /* peak memory usage for both UMFPACK_*symbolic and UMFPACK_numeric. */
    Info [UMFPACK_PEAK_MEMORY + what] =
	MAX (Symbolic->peak_sym_usage, num_usage) ;

    Info [UMFPACK_FLOPS + what] = flops ;
    Info [UMFPACK_LNZ + what] = lnz ;
    Info [UMFPACK_UNZ + what] = unz ;
}
