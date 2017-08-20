/* ========================================================================== */
/* === UMF_kernel_init_usage ================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* Return number of Units needed for initial elements by UMF_kernel_init */
/* (including the tuples). */

#include "umf_internal.h"
#include "umf_build_tuples_usage.h"

GLOBAL Int UMF_kernel_init_usage
(
    const Int Ap [ ],
    const Int Row_degree [ ],
    Int n_row,
    Int n_col,
    double *dusage		/* input and output argument */
)
{
    Int col, cdeg, usage ;

    /* elements: */
    usage = 0 ;
    for (col = 0 ; col < n_col ; col++)
    {
	DEBUG2 (("col "ID" Ap["ID"] = "ID"  Ap ["ID"] "ID"\n",
	    col, col+1, Ap [col+1], col, Ap [col])) ;
	cdeg = Ap [col+1] - Ap [col] ;
	ASSERT (cdeg >= 0) ;
	if (cdeg > 0)
	{
	    usage   += GET_ELEMENT_SIZE (cdeg, 1) + 1 ;
	    *dusage += DGET_ELEMENT_SIZE (cdeg, 1) + 1 ;
	}
    }

    DEBUG0 (("UMF_kernel_init_usage : original elements: "ID" %g\n",
	usage, *dusage)) ;

    /* tuples: */
    usage += UMF_build_tuples_usage ((Int *) NULL, (Int *) NULL,
	Row_degree, Row_degree, n_row, n_col, dusage) ;

    DEBUG0 (("UMF_kernel_init_usage : all: "ID" %g\n", usage, *dusage)) ;
    return (usage) ;
}

