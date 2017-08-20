/* ========================================================================== */
/* === umfpack_qsymbolic ==================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_qsymbolic
(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    const int Qinit [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

long umfpack_dl_qsymbolic
(
    long n_row,
    long n_col,
    const long Ap [ ],
    const long Ai [ ],
    const long Qinit [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

int umfpack_zi_qsymbolic
(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    const int Qinit [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

long umfpack_zl_qsymbolic
(
    long n_row,
    long n_col,
    const long Ap [ ],
    const long Ai [ ],
    const long Qinit [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    int n_row, n_col, *Ap, *Ai, *Qinit, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    status = umfpack_di_qsymbolic (n_row, n_col, Ap, Ai, Qinit, &Symbolic,
	Control, Info) ;

double long Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    long n_row, n_col, *Ap, *Ai, *Qinit, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    status = umfpack_dl_qsymbolic (n_row, n_col, Ap, Ai, Qinit, &Symbolic,
	Control, Info) ;

complex int Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    int n_row, n_col, *Ap, *Ai, *Qinit, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    status = umfpack_zi_qsymbolic (n_row, n_col, Ap, Ai, Qinit, &Symbolic,
	Control, Info) ;

complex long Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    long n_row, n_col, *Ap, *Ai, *Qinit, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    status = umfpack_zl_qsymbolic (n_row, n_col, Ap, Ai, Qinit, &Symbolic,
	Control, Info) ;

Purpose:

    Given the nonzero pattern of a sparse matrix A in column-oriented form, and
    a sparsity preserving column preordering Qinit, umfpack_*_qsymbolic performs
    the symbolic factorization of A*Qinit (or A (:,Qinit) in MATLAB notation).
    It also computes the column elimination tree post-ordering.  This is
    identical to umfpack_*_symbolic, except that colamd is not called and the
    user input column order Qinit is used instead.  Note that in general, the
    Qinit passed to umfpack_*_qsymbolic will differ from the final Q found in
    umfpack_*_numeric, because of the column etree postordering done in
    umfpack_*_qsymbolic and sparsity-preserving modifications made within each
    frontal matrix during umfpack_*_numeric.

    *** WARNING ***  A poor choice of Qinit can easily cause umfpack_*_numeric
    to use a huge amount of memory and do a lot of work.  The "default" symbolic
    analysis method is umfpack_*_symbolic, not this routine.  If you use this
    routine, the performance of UMFPACK is your responsibility;  UMFPACK will
    not try to second-guess a poor choice of Qinit.  If you are unsure about
    the quality of your Qinit, then call both umfpack_*_symbolic and
    umfpack_*_qsymbolic, and pick the one with lower estimates of work and
    memory usage (Info [UMFPACK_FLOPS_ESTIMATE] and
    Info [UMFPACK_PEAK_MEMORY_ESTIMATE]).  Don't forget to call
    umfpack_*_free_symbolic to free the Symbolic object that you don't need.

Returns:

    The value of Info [UMFPACK_STATUS]; see below.

Arguments:

    All arguments are the same as umfpack_*_symbolic, except for the following:

    Int Qinit [n_col] ;		Input argument, not modified.

	The user's fill-reducing initial column preordering.  This must be a
	permutation of 0..n_col-1.  If Qinit [k] = j, then column j is the kth
	column of the matrix A (:,Qinit) to be factorized.  If Qinit is an
	(Int *) NULL pointer, then colamd is called instead.  In fact,

	Symbolic = umfpack_*_symbolic (n_row, n_col, Ap, Ai, Control, Info) ;

	is identical to

	Symbolic = umfpack_*_qsymbolic (n_row, n_col, Ap, Ai, (Int *) NULL,
	    Control, Info) ;

    double Control [UMFPACK_CONTROL] ;	Input argument, not modified.

	Identical to umfpack_*_symbolic if Qinit is (Int *) NULL.  Otherwise,
	if Qinit is present, it is identical to umfpack_*_symbolic except for
	the following:

	Control [UMFPACK_DENSE_ROW]:  ignored.

	Control [UMFPACK_DENSE_COL]:  ignored.

    double Info [UMFPACK_INFO] ;	Output argument, not defined on input.

	Identical to umfpack_*_symbolic if Qinit is (Int *) NULL.  Otherwise,
	if Qinit is present, it is identical to umfpack_*_symbolic except for
	the following:

	Info [UMFPACK_NDENSE_ROW]:  zero
	Info [UMFPACK_NEMPTY_ROW]:  number of empty rows.
	Info [UMFPACK_NDENSE_COL]:  zero
	Info [UMFPACK_NEMPTY_COL]:  number of empty columns.
*/
