/* ========================================================================== */
/* === UMFPACK_symbolic ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Performs a symbolic factorization.
    See umfpack_symbolic.h for details.
*/

#include "umf_internal.h"

GLOBAL Int UMFPACK_symbolic
(
    Int n_row,
    Int n_col,
    const Int Ap [ ],
    const Int Ai [ ],
    void **SymbolicHandle,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
)
{
    Int *Qinit = (Int *) NULL ;
    return (UMFPACK_qsymbolic (n_row, n_col, Ap, Ai, Qinit, SymbolicHandle,
	Control, Info)) ;
}

