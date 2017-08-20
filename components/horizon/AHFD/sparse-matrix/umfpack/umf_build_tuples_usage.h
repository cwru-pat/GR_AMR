/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

GLOBAL Int UMF_build_tuples_usage
(
    const Int Col_tlen [ ],
    const Int Col_degree [ ],
    const Int Row_tlen [ ],
    const Int Row_degree [ ],
    Int n_row,
    Int n_col,
    double *dusage
) ;

