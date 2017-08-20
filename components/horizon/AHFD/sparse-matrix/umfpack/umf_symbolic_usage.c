/* ========================================================================== */
/* === UMF_symbolic_usage =================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* Returns the final size of the Symbolic object, in Units */

#include "umf_internal.h"

GLOBAL double UMF_symbolic_usage
(
    Int n_row,
    Int n_col,
    Int nchains,
    Int nfr
)
{
    double units =
	DUNITS (SymbolicType, 1)	/* Symbolic structure */
	+ DUNITS (Int, n_col+1)		/* Cperm_init */
	+ DUNITS (Int, n_row+1)		/* Rperm_init */
	+ 3*DUNITS (Int, nchains+1)	/* Chain_ */
	+ 4*DUNITS (Int, nfr+1) ;	/* Front_ */

    return (units) ;
}

