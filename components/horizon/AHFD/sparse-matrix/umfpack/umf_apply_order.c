/* ========================================================================== */
/* === UMF_apply_order ====================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Apply post-ordering of supernodal elimination tree.
*/

#include "umf_internal.h"

GLOBAL void UMF_apply_order
(
    Int Front [ ],
    const Int Order [ ],
    Int Temp [ ],
    Int n_col,
    Int nfr
)
{
    Int i, k ;
    for (i = 0 ; i < n_col ; i++)
    {
	k = Order [i] ;
	if (k != EMPTY)
	{
	    Temp [k] = Front [i] ;
	}
    }

    for (k = 0 ; k < nfr ; k++)
    {
	Front [k] = Temp [k] ;
    }
}

