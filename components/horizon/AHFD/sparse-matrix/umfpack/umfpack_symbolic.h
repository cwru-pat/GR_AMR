/* ========================================================================== */
/* === umfpack_symbolic ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_symbolic
(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

long umfpack_dl_symbolic
(
    long n_row,
    long n_col,
    const long Ap [ ],
    const long Ai [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

int umfpack_zi_symbolic
(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

long umfpack_zl_symbolic
(
    long n_row,
    long n_col,
    const long Ap [ ],
    const long Ai [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    int n_row, n_col, *Ap, *Ai, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    status = umfpack_di_symbolic (n_row, n_col, Ap, Ai, &Symbolic, Control,
	Info) ;

double long Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    long n_row, n_col, *Ap, *Ai, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    status = umfpack_dl_symbolic (n_row, n_col, Ap, Ai, &Symbolic, Control,
	Info) ;

complex int Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    int n_row, n_col, *Ap, *Ai, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    status = umfpack_zi_symbolic (n_row, n_col, Ap, Ai, &Symbolic, Control,
	Info) ;

complex long Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    long n_row, n_col, *Ap, *Ai, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    status = umfpack_zl_symbolic (n_row, n_col, Ap, Ai, &Symbolic, Control,
	Info) ;

Purpose:

    Given nonzero pattern of a sparse matrix A in column-oriented form,
    umfpack_*_symbolic performs a column pre-ordering to reduce fill-in
    (using UMF_colamd, modified from colamd V2.0 for UMFPACK), and a symbolic
    factorization.  This is required before the matrix can be numerically
    factorized with umfpack_*_numeric.  If you wish to bypass the UMF_colamd
    pre-ordering and provide your own ordering, use umfpack_*_qsymbolic instead.

Returns:

    The status code is returned.  See Info [UMFPACK_STATUS], below.

Arguments:

    Int n_row ;		Input argument, not modified.
    Int n_col ;		Input argument, not modified.

	A is an n_row-by-n_col matrix.  Restriction: n_row > 0 and n_col > 0.

    Int Ap [n_col+1] ;	Input argument, not modified.

	Ap is an integer array of size n_col+1.  On input, it holds the
	"pointers" for the column form of the sparse matrix A.  Column j of
	the matrix A is held in Ai [(Ap [j]) ... (Ap [j+1]-1)].  The first
	entry, Ap [0], must be zero, and Ap [j] <= Ap [j+1] must hold for all
	j in the range 0 to n_col-1.  The value nz = Ap [n_col] is thus the
	total number of entries in the pattern of the matrix A.  nz must be
	greater than or equal to zero.

    Int Ai [nz] ;	Input argument, not modified, of size nz = Ap [n_col].

	The nonzero pattern (row indices) for column j is stored in
	Ai [(Ap [j]) ... (Ap [j+1]-1)].  The row indices in a given column j
	must be in ascending order, and no duplicate row indices may be present.
	Row indices must be in the range 0 to n_row-1 (the matrix is 0-based).
	See umfpack_*_triplet_to_col for how to sort the columns of a matrix
	and sum up the duplicate entries.  See umfpack_*_report_matrix for how
	to print the matrix A.

    void **Symbolic ;	Output argument.

	**Symbolic is the address of a (void *) pointer variable in the user's
	calling routine (see Syntax, above).  On input, the contents of this
	variable are not defined.  On output, this variable holds a (void *)
	pointer to the Symbolic object (if successful), or (void *) NULL if
	a failure occurred.

    double Control [UMFPACK_CONTROL] ;	Input argument, not modified.

	If a (double *) NULL pointer is passed, then the default control
	settings are used.  Otherwise, the settings are determined from the
	Control array.  See umfpack_*_defaults on how to fill the Control
	array with the default settings.  If Control contains NaN's, the
	defaults are used.  The following Control parameters are used:

	Control [UMFPACK_DENSE_COL]:  columns with more than
	    max (16, Control [UMFPACK_DENSE_COL] * 16 * sqrt (n_row))
	    entries are placed placed last in the column pre-ordering by
	    UMF_colamd.  Default: 0.2.

	Control [UMFPACK_DENSE_ROW]:  rows with more than
	    max (16, Control [UMFPACK_DENSE_ROW] * 16 * sqrt (n_col))
	    entries (after "dense" columns are removed) are ignored in the
	    column pre-ordering, UMF_colamd.  Default: 0.2.

	Control [UMFPACK_BLOCK_SIZE]:  the block size to use for Level-3 BLAS
	    in the subsequent numerical factorization (umfpack_*_numeric).
	    A value less than 1 is treated as 1.  Default: 24.  Modifying this
	    parameter affects when updates are applied to the working frontal
	    matrix, and can indirectly affect fill-in and operation count.
	    As long as the block size is large enough (8 or so), this parameter
	    has modest effect on performance.  In Version 3.0, this parameter
	    was an input to umfpack_*_numeric, and had a default value of 16.
	    On a Sun UltraSparc, a block size of 24 is better for larger
	    matrices (16 is better for smaller ones, but not by much).  In the
	    current version, it is required in the symbolic analysis phase, and
	    is thus an input to this phase.

    double Info [UMFPACK_INFO] ;	Output argument, not defined on input.

	Contains statistics about the symbolic analysis.  If a (double *) NULL
	pointer is passed, then no statistics are returned in Info (this is not
	an error condition).  The entire Info array is cleared (all entries set
	to -1) and then the following statistics are computed:

	Info [UMFPACK_STATUS]: status code.  This is also the return value,
	    whether or not Info is present.

	    UMFPACK_OK

		Each column of the input matrix contained row indices
		in increasing order, with no duplicates.  Only in this case
		does umfpack_*_symbolic compute a valid symbolic factorization.
		For the other cases below, no Symbolic object is created
		(*Symbolic is (void *) NULL).

	    UMFPACK_ERROR_jumbled_matrix

		Columns of input matrix were jumbled (unsorted columns or
		duplicate entries).

	    UMFPACK_ERROR_n_nonpositive

		n is less than or equal to zero.

	    UMFPACK_ERROR_nz_negative

		Number of entries in the matrix is negative.

	    UMFPACK_ERROR_Ap0_nonzero

		Ap [0] is nonzero.

	    UMFPACK_ERROR_col_length_negative

		A column has a negative number of entries.

	    UMFPACK_ERROR_row_index_out_of_bounds

		A row index is out of bounds.

	    UMFPACK_ERROR_out_of_memory

		Insufficient memory to perform the symbolic analysis.

	    UMFPACK_ERROR_argument_missing

		One or more required arguments is missing.

	    UMFPACK_ERROR_problem_too_large

		Problem is too large; memory usage estimate causes an integer
		overlow.  If you are using umfpack_*i_symbolic, try using
		the long versions instead, umfpack_*l_symbolic.

	    UMFPACK_ERROR_internal_error

		Something very serious went wrong.  This is a bug.
		Please contact the author (davis@cise.ufl.edu).

	Info [UMFPACK_NROW]:  the value of the input argument n_row.

	Info [UMFPACK_NCOL]:  the value of the input argument n_col.

	Info [UMFPACK_NZ]:  the number of entries in the input matrix
	    (Ap [n_col]).

	Info [UMFPACK_SIZE_OF_UNIT]:  the number of bytes in a Unit,
	    for memory usage statistics below.

	Info [UMFPACK_SIZE_OF_INT]:  the number of bytes in an int.

	Info [UMFPACK_SIZE_OF_LONG]:  the number of bytes in a long.

	Info [UMFPACK_SIZE_OF_POINTER]:  the number of bytes in a void *
	    pointer.

	Info [UMFPACK_SIZE_OF_ENTRY]:  the number of bytes in a numerical entry.

	Info [UMFPACK_NDENSE_ROW]:  number of "dense" rows in A.  These rows are
	    ignored when the column pre-ordering is computed in UMF_colamd.
	    If > 0, then the matrix had to be re-analyzed by UMF_analyze, which
	    does not ignore these rows.

	Info [UMFPACK_NEMPTY_ROW]:  number of "empty" rows in A.  These are
	    rows that either have no entries, or whose entries are all in
	    "dense" columns.  Any given row is classified as either "dense"
	    or "empty" or "sparse".

	Info [UMFPACK_NDENSE_COL]:  number of "dense" columns in A.  These
	    columns are ordered last in the factorization, but before "empty"
	    columns.  Any given column is classified as either "dense" or
	    "empty" or "sparse".

	Info [UMFPACK_NEMPTY_COL]:  number of "empty" columns in A.  These are
	    columns that either have no entries, or whose entries are all in
	    "dense" rows.  These columns are ordered last in the factorization,
	    to the right of "dense" columns.

	Info [UMFPACK_SYMBOLIC_DEFRAG]:  number of garbage collections
	    performed in UMF_colamd, the column pre-ordering routine, and in
	    UMF_analyze, which is called if UMF_colamd isn't, or if UMF_colamd
	    ignores one or more "dense" rows.

	Info [UMFPACK_SYMBOLIC_PEAK_MEMORY]:  the amount of memory (in Units)
	    required for umfpack_*_symbolic to complete.  This is roughly
	    2.2*nz + (26 to 31)*n integers for a square matrix, depending on the
	    matrix.  This count includes the size of the Symbolic object itself,
	    which is reported in Info [UMFPACK_SYMBOLIC_SIZE].

	Info [UMFPACK_SYMBOLIC_SIZE]: the final size of the Symbolic object (in
	    Units).  This is fairly small, roughly 2*n to 9*n integers,
	    depending on the matrix.

	Info [UMFPACK_VARIABLE_INIT_ESTIMATE]: the Numeric object contains two
	    parts.  The first is fixed in size (O (n_row+n_col)).  The
	    second part holds the sparse LU factors and the contribution blocks
	    from factorized frontal matrices.  This part changes in size during
	    factorization.  Info [UMFPACK_VARIABLE_INIT_ESTIMATE] is the exact
	    size (in Units) required for this second variable-sized part in
	    order for the numerical factorization to start.

	Info [UMFPACK_VARIABLE_PEAK_ESTIMATE]: the estimated peak size (in
	    Units) of the variable-sized part of the Numeric object.  This is
	    usually an upper bound, but that is not guaranteed.

	Info [UMFPACK_VARIABLE_FINAL_ESTIMATE]: the estimated final size (in
	    Units) of the variable-sized part of the Numeric object.  This is
	    usually an upper bound, but that is not guaranteed.  It holds just
	    the sparse LU factors.

	Info [UMFPACK_NUMERIC_SIZE_ESTIMATE]:  an estimate of the final size (in
	    Units) of the entire Numeric object (both fixed-size and variable-
	    sized parts), which holds the LU factorization (including the L, U,
	    P and Q matrices).

	Info [UMFPACK_PEAK_MEMORY_ESTIMATE]:  an estimate of the total amount of
	    memory (in Units) required by umfpack_*_symbolic and
	    umfpack_*_numeric to perform both the symbolic and numeric
	    factorization.  This is the larger of the amount of memory needed
	    in umfpack_*_numeric itself, and the amount of memory needed in
	    umfpack_*_symbolic (Info [UMFPACK_SYMBOLIC_PEAK_MEMORY]).  The count
	    includes the size of both the Symbolic and Numeric objects
	    themselves.

	Info [UMFPACK_FLOPS_ESTIMATE]:  an estimate of the total floating-point
	    operations required to factorize the matrix.  This is a "true"
	    theoretical estimate of the number of flops that would be performed
	    by a flop-parsimonious sparse LU algorithm.  It assumes that no
	    extra flops are performed except for what is strictly required to
	    compute the LU factorization.  It ignores, for example, the flops
	    performed by umfpack_*_numeric to add contribution blocks of frontal
	    matrices together.  If L and U are the upper bound on the pattern
	    of the factors, then this flop count estimate can be represented in
	    MATLAB (for real matrices, not complex) as:

		Lnz = full (sum (spones (L))) - 1 ;	% nz in each col of L
		Unz = full (sum (spones (U')))' - 1 ;	% nz in each row of U
		flops = 2*Lnz*Unz + sum (Lnz) ;

	    The actual "true flop" count found by umfpack_*_numeric will be less
	    than this estimate.

	    For the real version, only (+ - * /) are counted.  For the complex
	    version, the following counts are used:

		operation	flops
	    	c = 1/b		6
		c = a*b		6
		c -= a*b	8

	Info [UMFPACK_LNZ_ESTIMATE]:  an estimate of the number of nonzeros in
	    L, including the diagonal.  Since L is unit-diagonal, the diagonal
	    of L is not stored.  This estimate is a strict upper bound on the
	    actual nonzeros in L to be computed by umfpack_*_numeric.

	Info [UMFPACK_UNZ_ESTIMATE]:  an estimate of the number of nonzeros in
	    U, including the diagonal.  This estimate is a strict upper bound on
	    the actual nonzeros in U to be computed by umfpack_*_numeric.

	Info [UMFPACK_SYMBOLIC_TIME]:  The time taken by umfpack_*_symbolic, in
	    seconds.  In the ANSI C version, this may be invalid if the time
	    taken is more than about 36 minutes, because of wrap-around in the
	    ANSI C clock function.  Compile UMFPACK with -DGETRUSAGE if you have
	    the more accurate getrusage function.

	At the start of umfpack_*_symbolic, all of Info is set of -1, and then
	after that only the above listed Info [...] entries are accessed.
	Future versions might modify different parts of Info.
*/
