/* ========================================================================== */
/* === UMF_build_tuples_usage =============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* Return number of Units needed for UMF_build_tuples */

#include "umf_internal.h"

GLOBAL Int UMF_build_tuples_usage
(
    const Int Col_tlen [ ],
    const Int Col_degree [ ],
    const Int Row_tlen [ ],
    const Int Row_degree [ ],
    Int n_row,
    Int n_col,
    double *dusage		/* input and output argument */
)
{
    Int row, col, usage ;
    double du ;

    /* note: tuple lengths are initialized, but the tuple lists themselves */
    /* may not be. */

    usage = 0 ;
    du = 0 ;
    if (!Col_tlen || !Col_degree)
    {
	/* Col_tlen and Col_degree arrays are missing, so this is the */
	/* initial matrix, with one element per column. */
	usage += n_col * (1 + UNITS (Tuple, 4)) ;
	du += ((double) n_col) * (1 + DUNITS (Tuple, 4)) ;
    }
    else
    {
	for (col = 0 ; col < n_col ; col++)
	{
	    if (NON_PIVOTAL_COL (col))
	    {
		usage += 1 + UNITS (Tuple, MAX (4, Col_tlen [col] + 1)) ;
		du += 1 + DUNITS (Tuple, MAX (4, Col_tlen [col] + 1)) ;
	    }
	}
    }
    ASSERT (Row_tlen && Row_degree) ;
    for (row = 0 ; row < n_row ; row++)
    {
	if (NON_PIVOTAL_ROW (row))
	{
	    usage += 1 + UNITS (Tuple, MAX (4, Row_tlen [row] + 1)) ;
	    du += 1 + DUNITS (Tuple, MAX (4, Row_tlen [row] + 1)) ;
	}
    }

    /* roundoff error in du is at most 3*n*epsilon */
    /* (here, and in UMF_kernel_init_usage) */
    /* Ignore NaN or Inf behavior here. */
    du = MAX (du, (double) usage * (1.0 + MAX_EPSILON)) ;
    du += (((double) n_row) + 2*((double) n_col)) * MAX_EPSILON ;
    du = ceil (du) ;

    DEBUG0 (("UMF_build_tuples_usage "ID" %g\n", usage, *dusage)) ;

    *dusage += du ;
    return (usage) ;
}

