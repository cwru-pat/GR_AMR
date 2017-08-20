/* ========================================================================== */
/* === UMFPACK_qsymbolic ==================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Performs a symbolic factorization.
    See umfpack_qsymbolic.h and umfpack_symbolic.h for details.

    Dynamic memory usage of UMFPACK_qsymbolic:

    1)  calls UMF_malloc 6 times, for workspace of total size (C + 5*n_col)
	integers.  The value of C is determined by the macro
	UMF_COLAMD_RECOMMENDED.  It is roughly
	max (2*nz, 3*n_col) + 8*n_col + 6*n_row + n_col + nz/5, or typically
	about 2.2*nz + 9*n_col + 6*n_row, where nz is the number of entries in
	A.  If A is square, this is 2.2*nz + 15*n.

    2)  three more calls to UMF_malloc are made, for a total space of
	(n_row+1 + n_col+1) integers + sizeof (SymbolicType).  
	sizeof (SymbolicType) is a small constant.  This space is part of the
	Symbolic object and is not freed unless an error occurs.  The analysis
	(UMF_colamd and/or UMF_analyze) is then performed.
	This is about 2n if A is square.

    3)  UMF_malloc is called 7 times, for a total workspace of
	(4*(nfr+1)+3*(nchains+1)) integers, where nfr is the total number of
	frontal matrices and nchains is the total number of frontal matrix
	chains, and where nchains <= nfr <= n_col.  This space is part of the
	Symbolic object and is not free'd unless an error occurs.  This is
	between 7 and about 7n integers when A is square.

	Thus, the peak memory usage of UMFPACK_symbolic occurs at this point.
	For a square matrix, it is at most about (2.2*nz + 24*n) integers for
	temporary workspace, plus between 2*n and 9*n integers for the final
	Symbolic object.

    4)	The 6 workspace objects allocated in step (1) are free'd via
	UMF_free.  The final Symbolic object consists of 10 allocated objects.
	Its final total size is lies roughly between 2*n and 9*n for a square
	matrix.  If an error occurs, all 10 objects are free'd via UMF_free.

    Dynamic memory usage of UMFPACK_free_symbolic:

    1)  All 10 objects comprising the Symbolic object are free'd via UMF_free.

*/

/* ========================================================================== */

#include "umf_internal.h"
#include "umf_symbolic_usage.h"
#include "umf_colamd.h"
#include "umf_set_stats.h"
#include "umf_analyze.h"
#include "umf_transpose.h"
#include "umf_is_permutation.h"
#include "umf_kernel_init_usage.h"
#include "umf_malloc.h"
#include "umf_free.h"

typedef struct	/* SWorkType */
{
    Int *Front_npivcol ;
    Int *Front_nrows ;
    Int *Front_ncols ;
    Int *Front_parent ;
    Int *Front_cols ;		/* in UMF_colamd only */
    Int *Ci ;

} SWorkType ;

PRIVATE void free_work
(
    SWorkType *SWork
) ;

PRIVATE void error
(
    SymbolicType **Symbolic,
    SWorkType *SWork
) ;

#define SYM_WORK_USAGE(n_col,Clen) (DUNITS (Int, Clen) + 5*DUNITS (Int, n_col))

/* required size of Ci for code that calls UMF_transpose and UMF_analyze below*/
#define UMF_ANALYZE_CLEN(nz,n_row,n_col,nn) \
    ((n_col) + MAX ((nz),1) + 3*(nn)+1 + (n_col))

#ifndef NDEBUG
PRIVATE Int init_count ;
#endif

/* ========================================================================== */

GLOBAL Int UMFPACK_qsymbolic
(
    Int n_row,
    Int n_col,
    const Int Ap [ ],
    const Int Ai [ ],
    const Int Q [ ],
    void **SymbolicHandle,
    const double Control [UMFPACK_CONTROL],
    double User_Info [UMFPACK_INFO]
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int colamd_ok, i, nz, maxfrsize, s, j, newj, status, maxoldpiv, chain_npiv,
	maxnrows, maxncols, nfr, col, nchains, maxrows, maxcols, p, nb, nn,
	*Chain_start, *Chain_maxrows, *Chain_maxcols, *Front_npivcol, *Ci,
	Clen, colamd_stats [COLAMD_STATS], fpiv, frows, fcols, n_inner,
	child, cr, cc, cp, parent, *Link, *Row_degree, row, *Front_parent,
	head_usage, tail_usage, max_usage, analyze_compactions,
	lf, uf, init_tail_usage, k, do_colamd, do_UMF_analyze, too_large,
	fpivcol, fallrows, fallcols, *InFront, *F1,
	*Front_1strow, f1rows, kk, *Cperm_init, *Rperm_init, newrow,
	*Front_leftmostdesc, Clen_analyze ;
    double knobs [COLAMD_KNOBS], flops, f, r, c, tstart, tend,
	*Info, Info2 [UMFPACK_INFO], drow, dcol, dusage, dlf, duf, dmax_usage,
	dhead_usage, dh, dlnz, dunz, dmaxfrsize, dClen, dClen_analyze ;
    Int nempty_row, nempty_col ;
    SymbolicType *Symbolic ;
    SWorkType SWorkSpace, *SWork ;

#ifndef NDEBUG
    double f2 ;
    Int lf2, uf2 ;
    UMF_dump_start ( ) ;
    init_count = UMF_malloc_count ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get the amount of time used by the process so far */
    /* ---------------------------------------------------------------------- */

    tstart = umfpack_timer ( ) ;

    /* ---------------------------------------------------------------------- */
    /* check input parameters */
    /* ---------------------------------------------------------------------- */

    if (Control)
    {
	/* use the Control array passed to us by the caller. */

	if (SCALAR_IS_NAN (Control [UMFPACK_DENSE_ROW]))
	{
	    drow = UMFPACK_DEFAULT_DENSE_ROW ;
	}
	else
	{
	    drow = Control [UMFPACK_DENSE_ROW] ;
	}

	if (SCALAR_IS_NAN (Control [UMFPACK_DENSE_COL]))
	{
	    dcol = UMFPACK_DEFAULT_DENSE_COL ;
	}
	else
	{
	    dcol = Control [UMFPACK_DENSE_COL] ;
	}

	if (SCALAR_IS_NAN (Control [UMFPACK_BLOCK_SIZE]))
	{
	    nb = UMFPACK_DEFAULT_BLOCK_SIZE ;
	}
	else
	{
	    nb = (Int) Control [UMFPACK_BLOCK_SIZE] ;
	}

    }
    else
    {
	/* no Control passed - use defaults instead */
	drow = UMFPACK_DEFAULT_DENSE_ROW ;
	dcol = UMFPACK_DEFAULT_DENSE_COL ;
	nb = UMFPACK_DEFAULT_BLOCK_SIZE ;
    }

    nb = MAX (1, nb) ;
    DEBUG0 (("nb = "ID"\n", nb)) ;

    if (User_Info)
    {
	/* return Info in user's array */
	Info = User_Info ;
    }
    else
    {
	/* no Info array passed - use local one instead */
	Info = Info2 ;
    }
    for (i = 0 ; i < UMFPACK_INFO ; i++)
    {
	Info [i] = EMPTY ;
    }

    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;

    Info [UMFPACK_STATUS] = UMFPACK_OK ;
    Info [UMFPACK_NROW] = n_row ;
    Info [UMFPACK_NCOL] = n_col ;
    Info [UMFPACK_SIZE_OF_UNIT] = (double) (sizeof (Unit)) ;
    Info [UMFPACK_SIZE_OF_INT] = (double) (sizeof (int)) ;
    Info [UMFPACK_SIZE_OF_LONG] = (double) (sizeof (long)) ;
    Info [UMFPACK_SIZE_OF_POINTER] = (double) (sizeof (void *)) ;
    Info [UMFPACK_SIZE_OF_ENTRY] = (double) (sizeof (Entry)) ;
    Info [UMFPACK_SYMBOLIC_DEFRAG] = 0 ;

    if (!Ai || !Ap || !SymbolicHandle)
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_argument_missing ;
	return (UMFPACK_ERROR_argument_missing) ;
    }

    *SymbolicHandle = (void *) NULL ;

    if (n_row <= 0 || n_col <= 0)	/* n_row, n_col must be > 0 */
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_n_nonpositive ;
	return (UMFPACK_ERROR_n_nonpositive) ;
    }

    nz = Ap [n_col] ;
    DEBUG0 (("In UMFPACK_symbolic, n_row "ID" n_col "ID" nz "ID"\n",
	n_row, n_col, nz)) ;
    Info [UMFPACK_NZ] = nz ;
    if (nz < 0)
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_nz_negative ;
	return (UMFPACK_ERROR_nz_negative) ;
    }

    if (Q)
    {
	do_colamd = FALSE ;
    }
    else
    {
	do_colamd = TRUE ;
    }

    /* ---------------------------------------------------------------------- */
    /* determine amount of memory required for UMFPACK_symbolic */
    /* ---------------------------------------------------------------------- */

    /* The size of Clen required for UMF_colamd is always larger than */
    /* UMF_analyze, but the max is included here in case that changes in */
    /* future versions. */

    dClen = UMF_COLAMD_RECOMMENDED ((double) nz, (double) n_row,
	(double) n_col) ;
    dClen_analyze = UMF_ANALYZE_CLEN ((double) nz, (double) n_row,
	(double) n_col, (double) nn) ;
    dClen = MAX (dClen, dClen_analyze) ;
    too_large = INT_OVERFLOW (dClen * sizeof (Int)) ;

    if (too_large)
    {
	/* Problem is too large for array indexing (Ci [i]) with an Int i. */
	/* Cannot even analyze the problem to determine upper bounds on */
	/* memory usage. Need to use the long integer version, umfpack_*l_*. */
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_problem_too_large ;
	Info [UMFPACK_SYMBOLIC_PEAK_MEMORY] =
	    SYM_WORK_USAGE (n_col, dClen) +
	    UMF_symbolic_usage (n_row, n_col, n_col, n_col) ;
	return (UMFPACK_ERROR_problem_too_large) ;
    }

    Clen = UMF_COLAMD_RECOMMENDED (nz, n_row, n_col) ;
    Clen_analyze = UMF_ANALYZE_CLEN (nz, n_row, n_col, nn) ;
    Clen = MAX (Clen, Clen_analyze) ;

    /* worst case total memory usage for UMFPACK_symbolic (revised below) */
    Info [UMFPACK_SYMBOLIC_PEAK_MEMORY] =
	SYM_WORK_USAGE (n_col, Clen) +
	UMF_symbolic_usage (n_row, n_col, n_col, n_col) ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    Symbolic = (SymbolicType *) NULL ;
    SWork = &SWorkSpace ;	/* used for UMFPACK_symbolic only */

    /* Note that SWork->Front_* does not include the dummy placeholder front. */
    /* This space is accounted for by the SYM_WORK_USAGE (n_col, Clen) macro. */
    SWork->Ci = (Int *) UMF_malloc (Clen, sizeof (Int)) ;
    SWork->Front_npivcol = (Int *) UMF_malloc (n_col, sizeof (Int)) ;
    SWork->Front_nrows = (Int *) UMF_malloc (n_col, sizeof (Int)) ;
    SWork->Front_ncols = (Int *) UMF_malloc (n_col, sizeof (Int)) ;
    SWork->Front_parent = (Int *) UMF_malloc (n_col, sizeof (Int)) ;
    SWork->Front_cols = (Int *) UMF_malloc (n_col, sizeof (Int)) ;

    if (!SWork->Ci || !SWork->Front_npivcol || !SWork->Front_nrows ||
	!SWork->Front_ncols || !SWork->Front_parent || !SWork->Front_cols)
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	error (&Symbolic, SWork) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate the first part of the Symbolic object (header and Cperm_init) */
    /* ---------------------------------------------------------------------- */

    Symbolic = (SymbolicType *) UMF_malloc (1, sizeof (SymbolicType)) ;

    if (!Symbolic)
    {
	/* If we fail here, Symbolic is NULL and thus it won't be */
	/* dereferenced by UMFPACK_free_symbolic, as called by error ( ). */
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	error (&Symbolic, SWork) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }

    /* We now know that Symbolic has been allocated */
    Symbolic->valid = 0 ;
    Symbolic->Chain_start = (Int *) NULL ;
    Symbolic->Chain_maxrows = (Int *) NULL ;
    Symbolic->Chain_maxcols = (Int *) NULL ;
    Symbolic->Front_npivcol = (Int *) NULL ;
    Symbolic->Front_parent = (Int *) NULL ;
    Symbolic->Front_1strow = (Int *) NULL ;
    Symbolic->Front_leftmostdesc = (Int *) NULL ;

    Symbolic->Cperm_init = (Int *) UMF_malloc (n_col+1, sizeof (Int)) ;
    Symbolic->Rperm_init = (Int *) UMF_malloc (n_row+1, sizeof (Int)) ;
    Cperm_init = Symbolic->Cperm_init ;
    Rperm_init = Symbolic->Rperm_init ;

    if (!Symbolic->Cperm_init || !Symbolic->Rperm_init)
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	error (&Symbolic, SWork) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }

    DEBUG0 (("(1)Symbolic UMF_malloc_count - init_count = "ID"\n",
	UMF_malloc_count - init_count)) ;
    ASSERT (UMF_malloc_count == init_count + 9) ;

    Symbolic->n_row = n_row ;
    Symbolic->n_col = n_col ;
    Symbolic->nz = nz ;
    Symbolic->drow = drow ;
    Symbolic->dcol = dcol ;
    Symbolic->nb = nb ;

    /* ---------------------------------------------------------------------- */
    /* use colamd or input column ordering Q */
    /* ---------------------------------------------------------------------- */

    if (do_colamd)
    {

	/* ------------------------------------------------------------------ */
	/* copy the matrix into colamd workspace (colamd destroys its input) */
	/* ------------------------------------------------------------------ */

	Ci = SWork->Ci ;
	for (i = 0 ; i < nz ; i++)
	{
	    DEBUG5 (("Ai  "ID": "ID"\n", i, Ai [i])) ;
	    Ci [i] = Ai [i] ;
	}

	/* load the column pointers into Cperm_init [0..n_col] */
	for (col = 0 ; col <= n_col ; col++)
	{
	    DEBUG5 (("Ap  "ID": "ID"\n", col, Ap [col])) ;
	    Cperm_init [col] = Ap [col] ;
	}

	/* ------------------------------------------------------------------ */
	/* set UMF_colamd defaults */
	/* ------------------------------------------------------------------ */

	UMF_colamd_set_defaults (knobs) ;
	knobs [COLAMD_DENSE_ROW] = drow ;
	knobs [COLAMD_DENSE_COL] = dcol ;

	/* ------------------------------------------------------------------ */
	/* check input matrix and find the initial column pre-ordering */
	/* ------------------------------------------------------------------ */

	colamd_ok = UMF_colamd (n_row, n_col, Clen, SWork->Ci,
		Cperm_init, knobs, colamd_stats,
		SWork->Front_npivcol,
		SWork->Front_nrows,
		SWork->Front_ncols,
		SWork->Front_parent,
		/* used in UMF_colamd only: */ SWork->Front_cols,
		&nfr) ;

	/* get number of truly empty rows and columns */
	nempty_row = colamd_stats [COLAMD_EMPTY_ROW] ;
	nempty_col = colamd_stats [COLAMD_EMPTY_COL] ;

	if (!colamd_ok)
	{
	    /* This will be modified below.  It is not (yet) an internal */
	    /* error.  It is only an internal error if is it not modified */
	    /* below (colamd_stats missing). */
	    Info [UMFPACK_STATUS] = UMFPACK_ERROR_internal_error ;
	}

	s = colamd_stats [COLAMD_STATUS] ;

	if (s == COLAMD_OK)
	{
	    status = UMFPACK_OK ;
	}
	else if (s == COLAMD_ERROR_jumbled_matrix)
	{
	    status = UMFPACK_ERROR_jumbled_matrix ;
	}
	else if (s == COLAMD_ERROR_p0_nonzero)
	{
	    status = UMFPACK_ERROR_Ap0_nonzero ;
	}
	else if (s == COLAMD_ERROR_row_index_out_of_bounds)
	{
	    status = UMFPACK_ERROR_row_index_out_of_bounds ;
	}
	else if (s == COLAMD_ERROR_col_length_negative)
	{
	    status = UMFPACK_ERROR_col_length_negative ;
	}
	else
	{
	    /* these errors "cannot" happen */
	    /* COLAMD_ERROR_A_too_small */
	    /* COLAMD_ERROR_A_not_present */
	    /* COLAMD_ERROR_p_not_present */
	    /* COLAMD_ERROR_out_of_memory */
	    /* COLAMD_ERROR_internal_error */
	    /* COLAMD_ERROR_nrow_negative */
	    /* COLAMD_ERROR_ncol_negative */
	    /* COLAMD_ERROR_nnz_negative */
	    status = UMFPACK_ERROR_internal_error ;
	}

	Info [UMFPACK_STATUS] = status ;
	if (status != UMFPACK_OK)
	{
	    error (&Symbolic, SWork) ;
	    return (status) ;
	}

	Info [UMFPACK_NDENSE_ROW] = colamd_stats [COLAMD_DENSE_ROW] ;
	Info [UMFPACK_NEMPTY_ROW] = colamd_stats [COLAMD_EMPTY_ROW]
		+ colamd_stats [COLAMD_NEWLY_EMPTY_ROW] ;
	Info [UMFPACK_NDENSE_COL] = colamd_stats [COLAMD_DENSE_COL] ;
	Info [UMFPACK_NEMPTY_COL] = colamd_stats [COLAMD_EMPTY_COL]
		+ colamd_stats [COLAMD_NEWLY_EMPTY_COL] ;
	Info [UMFPACK_SYMBOLIC_DEFRAG] = colamd_stats [COLAMD_DEFRAG_COUNT];

	/* re-analyze if any "dense" rows or cols ignored by UMF_colamd */
	do_UMF_analyze = 
	    colamd_stats [COLAMD_DENSE_ROW] > 0 ||
	    colamd_stats [COLAMD_DENSE_COL] > 0 ;

#ifndef NDEBUG
	for (col = 0 ; col < n_col ; col++)
	{
	    DEBUG1 (("Cperm_init ["ID"] = "ID"\n", col, Cperm_init[col]));
	}
	/* make sure colamd returned a valid permutation */
	ASSERT (Cperm_init) ;
	ASSERT (UMF_is_permutation (Cperm_init, SWork->Ci, n_col, n_col)) ;
#endif

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* do not call colamd - use input Q instead */
	/* ------------------------------------------------------------------ */

	Int length ;

	/* use Ci as workspace to check input permutation */
	if (!UMF_is_permutation (Q, SWork->Ci, n_col, n_col))
	{
	    Info [UMFPACK_STATUS] = UMFPACK_ERROR_invalid_permutation ;
	    error (&Symbolic, SWork) ;
	    return (UMFPACK_ERROR_invalid_permutation) ;
	}

	if (Ap [0] != 0)
	{
	    Info [UMFPACK_STATUS] = UMFPACK_ERROR_Ap0_nonzero ;
	    error (&Symbolic, SWork) ;
	    return (UMFPACK_ERROR_Ap0_nonzero) ;
	}

	nempty_col = 0 ;

	/* check Ap input */
	for (j = 0 ; j < n_col ; j++)
	{
	    length = Ap [j+1] - Ap [j] ;
	    if (length == 0)
	    {
	    	/* this is an empty column */
		DEBUG1 (("Original empty column: "ID"\n", j)) ;
		nempty_col++ ;
	    }
	    if (length < 0)
	    {
		Info [UMFPACK_STATUS] = UMFPACK_ERROR_col_length_negative ;
		error (&Symbolic, SWork) ;
		return (UMFPACK_ERROR_col_length_negative) ;
	    }
	}

	Info [UMFPACK_NDENSE_ROW] = 0 ;
	Info [UMFPACK_NEMPTY_ROW] = 0 ;	/* this is fixed below */
	Info [UMFPACK_NDENSE_COL] = 0 ;
	Info [UMFPACK_NEMPTY_COL] = nempty_col ;

	if (nempty_col == 0)
	{
	    /* copy the user's input permutation */
	    for (k = 0 ; k < n_col ; k++)
	    {
		Cperm_init [k] = Q [k] ;
	    }
	}
	else
	{
	    /* partition the user's input permutation */

	    /* move empty columns last */
	    k = n_col ;
	    for (j = n_col-1 ; j >= 0 ; j--)
	    {
		col = Q [j] ;
		length = Ap [col+1] - Ap [col] ;
		if (length == 0)
		{
		    Cperm_init [--k] = col ;
		    DEBUG1 (("Moving empty col "ID" last: "ID"\n", col, k)) ;
		}
	    }
	    ASSERT (n_col-k == nempty_col) ;

	    /* move remaining columns first */
	    k = 0 ;
	    for (j = 0 ; j < n_col ; j++)
	    {
		col = Q [j] ;
		length = Ap [col+1] - Ap [col] ;
		if (length != 0)
		{
		    DEBUG1 (("Non empty col "ID" at: "ID"\n", col, k)) ;
		    Cperm_init [k++] = col ;
		}
	    }
	}

	do_UMF_analyze = TRUE ;

    }

    Cperm_init [n_col] = EMPTY ;	/* unused in Cperm_init */

    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    DEBUG3 (("Cperm_init column permutation:\n")) ;
    ASSERT (UMF_is_permutation (Cperm_init, SWork->Ci, n_col, n_col)) ;
    for (k = 0 ; k < n_col ; k++)
    {
	DEBUG3 ((ID"\n", Cperm_init [k])) ;
    }
    /* ensure that empty columns have been placed last in A (:,Cperm_init) */
    for (newj = 0 ; newj < n_col ; newj++)
    {
	/* empty columns will be last in A (:, Cperm_init (1:n_col)) */
	j = Cperm_init [newj] ;
	ASSERT (IMPLIES (newj >= n_col-nempty_col, Ap [j+1] - Ap [j] == 0)) ;
	ASSERT (IMPLIES (newj <  n_col-nempty_col, Ap [j+1] - Ap [j] > 0)) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* analyze, using the given ordering, or using colamd's ordering */
    /* ---------------------------------------------------------------------- */

    if (do_UMF_analyze)
    {

	Int *W, *Bp, *Bi, *Cperm2, ok, ilast, *P, Clen2, bsize, Clen0 ;

	/* ------------------------------------------------------------------ */
	/* Ci [0 .. Clen-1] holds the following work arrays:

		first Clen0 entries	empty space, where Clen0 =
					Clen - (nn+1 + 2*nn + n_col)
					and Clen0 >= nz + n_col
		next nn+1 entries	Bp [0..nn]
		next nn entries		Link [0..nn-1]
		next nn entries		W [0..nn-1]
		last n_col entries	Cperm2 [0..n_col-1]
	*/

	Ci = SWork->Ci ;
	Clen0 = Clen - (nn+1 + 2*nn + n_col) ;
	Bp = Ci + Clen0 ;
	Link = Bp + (nn+1) ;
	W = Link + nn ;
	Cperm2 = W + nn ;
	ASSERT (Cperm2 + n_col == Ci + Clen) ;
	ASSERT (Clen0 >= nz + n_col) ;

	/* ------------------------------------------------------------------ */
	/* P = order that rows will be used in UMF_analyze */
	/* ------------------------------------------------------------------ */

	/* use W to mark rows, and use Link for row permutation P [ [ */
	for (row = 0 ; row < n_row ; row++)
	{
	    W [row] = FALSE ;
	}
	P = Link ;

	k = 0 ;
	for (newj = 0 ; newj < n_col ; newj++)
	{
	    /* empty columns will be last in A (:,Cperm_init) */
	    j = Cperm_init [newj] ;
	    ilast = -1 ;
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		row = Ai [p] ;
		if (row < 0 || row >= n_row)
		{
		    Info [UMFPACK_STATUS] =
			UMFPACK_ERROR_row_index_out_of_bounds ;
		    error (&Symbolic, SWork) ;
		    return (UMFPACK_ERROR_row_index_out_of_bounds) ;
		}
		if (row <= ilast)
		{
		    Info [UMFPACK_STATUS] = UMFPACK_ERROR_jumbled_matrix ;
		    error (&Symbolic, SWork) ;
		    return (UMFPACK_ERROR_jumbled_matrix) ;
		}
		if (!W [row])
		{
		    /* this row has just been see for the first time */
		    W [row] = TRUE ;
		    P [k++] = row ;
		}
		ilast = row ;
	    }
	}

	/* If the matrix has truly empty rows, then P will not be */
	/* complete, and visa versa.  The matrix is structurally singular. */
	ASSERT (IMPLIES (do_colamd, nempty_row == n_row-k)) ;
	nempty_row = n_row - k ;
	Info [UMFPACK_NEMPTY_ROW] = nempty_row ;
	if (k < n_row)
	{
	    /* complete P by putting empty rows last in their natural order, */
	    /* rather than declaring an error (the matrix is singular) */
	    for (row = 0 ; row < n_row ; row++)
	    {
		if (!W [row])
		{
		    /* W [row] = TRUE ;  (not required) */
		    P [k++] = row ;
		}
	    }
	}

	/* contents of W no longer needed ] */

#ifndef NDEBUG
	DEBUG3 (("Induced row permutation:\n")) ;
	ASSERT (k == n_row) ;
	ASSERT (UMF_is_permutation (P, W, n_row, n_row)) ;
	for (k = 0 ; k < n_row ; k++)
	{
	    DEBUG3 ((ID"\n", P [k])) ;
	}
#endif

	/* ------------------------------------------------------------------ */
	/* B = row-form of the pattern of A (P, Cperm_init (1:n_col)) */
	/* ------------------------------------------------------------------ */

	/* Ci [0 .. Clen-1] holds the following work arrays:

		first Clen2 entries	empty space, must be at least >= n_col
		next max (nz,1)		Bi [0..max (nz,1)-1]
		next nn+1 entries	Bp [0..nn]
		next nn entries		Link [0..nn-1]
		next nn entries		W [0..nn-1]
		last n_col entries	Cperm2 [0..n_col-1]

		This memory usage is accounted for by the UMF_ANALYZE_CLEN
		macro.
	*/

	Clen2 = Clen0 ;
	bsize = MAX (nz, 1) ;
	Clen2 -= bsize ;
	Bi = Ci + Clen2 ;
	ASSERT (Clen2 >= n_col) ;

	/* skip error test (already done, above).  Do pattern only. */
	(void) UMF_transpose (n_row, n_col, Ap, Ai, (double *) NULL, P,
	    Cperm_init, n_col - nempty_col, Bp, Bi, (double *) NULL, W, FALSE
#ifdef COMPLEX
	    , (double *) NULL, (double *) NULL, FALSE
#endif
	    ) ;

	/* contents of P (same as Link) and W not needed */
	/* still need Link and W as work arrays, though ] */

	ASSERT (Bp [0] == 0) ;
	ASSERT (Bp [n_row] == nz) ;

	/* increment Bp to point into Ci, not Bi */
	for (i = 0 ; i <= n_row ; i++)
	{
	    Bp [i] += Clen2 ;
	}
	ASSERT (Bp [0] == Clen0 - bsize) ;
	ASSERT (Bp [n_row] <= Clen0) ;

	/* Ci [0 .. Clen-1] holds the following work arrays:

		first Clen0 entries	Ci [0 .. Clen0-1], where the col indices
					of B are at the tail end of this part,
					and Bp [0] = Clen2 >= n_col.  Note
					that Clen0 = Clen2 + max (nz,1).
		next nn+1 entries	Bp [0..nn]
		next nn entries		Link [0..nn-1]
		next nn entries		W [0..nn-1]
		last n_col entries	Cperm2 [0..n_col-1]
	*/

	/* ------------------------------------------------------------------ */
	/* analyze */
	/* ------------------------------------------------------------------ */

	/* only analyze the non-empty part of the matrix */
	ok = UMF_analyze (n_row - nempty_row, n_col - nempty_col,
		Ci, Bp, Cperm2, W, Link,
		SWork->Front_ncols, SWork->Front_nrows, SWork->Front_npivcol,
		SWork->Front_parent, &nfr, &analyze_compactions) ;
	if (!ok)
	{
	    Info [UMFPACK_STATUS] = UMFPACK_ERROR_internal_error ;
	    error (&Symbolic, SWork) ;
	    return (UMFPACK_ERROR_internal_error) ;
	}
	Info [UMFPACK_SYMBOLIC_DEFRAG] += analyze_compactions ;

	/* ------------------------------------------------------------------ */
	/* combine the input permutation and UMF_analyze's permutation */
	/* ------------------------------------------------------------------ */

	/* Cperm2 is the column etree post-ordering */
	ASSERT (UMF_is_permutation (Cperm2, W,
	    n_col-nempty_col, n_col-nempty_col)) ;

	/* Note that the empty columns remain at the end of Cperm_init */
	for (k = 0 ; k < n_col - nempty_col ; k++)
	{
	    W [k] = Cperm_init [Cperm2 [k]] ;
	}

	for (k = 0 ; k < n_col - nempty_col ; k++)
	{
	    Cperm_init [k] = W [k] ;
	}

	ASSERT (UMF_is_permutation (Cperm_init, W, n_col, n_col)) ;

    }

    /* ---------------------------------------------------------------------- */
    /* determine the size of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    nchains = 0 ;
    for (i = 0 ; i < nfr ; i++)
    {
	if (SWork->Front_parent [i] != i+1)
	{
	    nchains++ ;
	}
    }

    Symbolic->nchains = nchains ;
    Symbolic->nfr = nfr ;

    /* final size of Symbolic object */
    Info [UMFPACK_SYMBOLIC_SIZE] =
	UMF_symbolic_usage (n_row, n_col, nchains, nfr) ;

    /* actual peak memory usage for UMFPACK_symbolic (actual nfr, nchains) */
    Info [UMFPACK_SYMBOLIC_PEAK_MEMORY] =
	SYM_WORK_USAGE (n_col, Clen) + Info [UMFPACK_SYMBOLIC_SIZE] ;
    Symbolic->peak_sym_usage = Info [UMFPACK_SYMBOLIC_PEAK_MEMORY] ;

    DEBUG0 (("Number of fronts: "ID"\n", nfr)) ;

    /* ---------------------------------------------------------------------- */
    /* allocate the second part of the Symbolic object (Front_*, Chain_*) */
    /* ---------------------------------------------------------------------- */

    /* Note that Symbolic->Front_* does include the dummy placeholder front */
    Symbolic->Front_npivcol = (Int *) UMF_malloc (nfr+1, sizeof (Int)) ;
    Symbolic->Front_parent = (Int *) UMF_malloc (nfr+1, sizeof (Int)) ;
    Symbolic->Front_1strow = (Int *) UMF_malloc (nfr+1, sizeof (Int)) ;
    Symbolic->Front_leftmostdesc = (Int *) UMF_malloc (nfr+1, sizeof (Int)) ;

    Symbolic->Chain_start = (Int *) UMF_malloc (nchains+1, sizeof (Int)) ;
    Symbolic->Chain_maxrows = (Int *) UMF_malloc (nchains+1, sizeof (Int)) ;
    Symbolic->Chain_maxcols = (Int *) UMF_malloc (nchains+1, sizeof (Int)) ;

    if (!Symbolic->Front_npivcol || !Symbolic->Front_parent ||
	!Symbolic->Front_1strow || !Symbolic->Front_leftmostdesc ||
	!Symbolic->Chain_start || !Symbolic->Chain_maxrows ||
	!Symbolic->Chain_maxcols)
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	error (&Symbolic, SWork) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }
    DEBUG0 (("(2)Symbolic UMF_malloc_count - init_count = "ID"\n",
	UMF_malloc_count - init_count)) ;
    ASSERT (UMF_malloc_count == init_count + 16) ;

    Front_npivcol = Symbolic->Front_npivcol ;
    Front_parent = Symbolic->Front_parent ;
    Front_1strow = Symbolic->Front_1strow ;
    Front_leftmostdesc = Symbolic->Front_leftmostdesc ;

    Chain_start = Symbolic->Chain_start ;
    Chain_maxrows = Symbolic->Chain_maxrows ;
    Chain_maxcols = Symbolic->Chain_maxcols ;

    /* ---------------------------------------------------------------------- */
    /* find row degrees and assign rows to fronts */
    /* ---------------------------------------------------------------------- */

    /* Use SWork->Ci as temporary workspace for Row_degree, InFront, and F1 */
    Row_degree = SWork->Ci ;		/* [ of size n_row */
    InFront = Row_degree + n_row ;	/* [ of size n_row */
    F1 = InFront + n_row ;		/* [ of size nfr+1 */
    ASSERT (Clen >= 2*n_row + nfr+1) ;
    for (row = 0 ; row < n_row ; row++)
    {
	Row_degree [row] = 0 ;
	InFront [row] = nfr ;	/* empty rows go to dummy front nfr */
    }

    newj = 0 ;
    k = 0 ;
    for (i = 0 ; i < nfr ; i++)
    {
	fpivcol = SWork->Front_npivcol [i] ;
	DEBUG1 (("Front "ID" k "ID" npivcol "ID" nrows "ID" ncols "ID"\n",
	    i, k, fpivcol, SWork->Front_nrows [i], SWork->Front_ncols [i])) ;
	k += fpivcol ;

	/* copy Front info into Symbolic object from SWork */
	Front_npivcol [i] = fpivcol ;
	Front_parent [i] = SWork->Front_parent [i] ;

	f1rows = 0 ;

	/* for all pivot columns in front i */
        for (kk = 0 ; kk < fpivcol ; kk++, newj++)
        {
	    j = Cperm_init [newj] ;
	    ASSERT (IMPLIES (newj >= n_col-nempty_col, Ap [j+1] - Ap [j] == 0));
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		row = Ai [p] ;
		if (Row_degree [row] == 0)
		{
		    /* this row belongs to front i */
		    DEBUG1 (("    Row "ID" in Front "ID"\n", row, i)) ;
		    InFront [row] = i ;
		    f1rows++ ;
		}
	        Row_degree [row]++ ;
	    }
    	}
	Front_1strow [i] = f1rows ;
    }

    /* assign empty columns to dummy placehold front nfr */
    DEBUG1 (("Dummy Cols in Front "ID" :: "ID"\n", nfr, n_col-k)) ;
    Front_npivcol [nfr] = n_col - k ;
    Front_parent [nfr] = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* find initial row permutation */
    /* ---------------------------------------------------------------------- */

    /* determine the first row in each front (in the new row ordering) */
    k = 0 ;
    for (i = 0 ; i < nfr ; i++)
    {
	f1rows = Front_1strow [i] ;
	DEBUG1 (("Front "ID" :: npivcol "ID" parent "ID,
	    i, Front_npivcol [i], Front_parent [i])) ;
	DEBUG1 (("    1st rows in Front "ID" : "ID"\n", i, f1rows)) ;
	Front_1strow [i] = k ;
    	k += f1rows ;
    }

    /* assign empty rows to dummy placehold front nfr */
    DEBUG1 (("Rows in Front "ID" (dummy): "ID"\n", nfr, n_row-k)) ;
    Front_1strow [nfr] = k ;
    DEBUG1 (("nfr "ID" 1strow[nfr] "ID" nrow "ID"\n", nfr, k, n_row)) ;

    for (i = 0 ; i <= nfr ; i++)
    {
	F1 [i] = Front_1strow [i] ;
    }

    for (row = 0 ; row < n_row ; row++)
    {
	i = InFront [row] ;
	newrow = F1 [i]++ ;
	Rperm_init [newrow] = row ;
    }
    Rperm_init [n_row] = EMPTY ;	/* unused */

#ifndef NDEBUG
    for (k = 0 ; k < n_row ; k++)
    {
    	DEBUG2 (("Rperm_init ["ID"] = "ID"\n", k, Rperm_init [k])) ;
    }
#endif

    /* ] done using F1 */
    /* ] done using InFront */

    /* ---------------------------------------------------------------------- */
    /* find the leftmost descendant of each front */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i <= nfr ; i++)
    {
	Front_leftmostdesc [i] = EMPTY ;
    }

    for (i = 0 ; i < nfr ; i++)
    {
	/* start at i and walk up the tree */
	DEBUG2 (("Walk up front tree from "ID"\n", i)) ;
	j = i ;
	while (j != EMPTY && Front_leftmostdesc [j] == EMPTY)
	{
	    DEBUG3 (("  Leftmost desc of "ID" is "ID"\n", j, i)) ;
	    Front_leftmostdesc [j] = i ;
	    j = Front_parent [j] ;
	    DEBUG3 (("  go to j = "ID"\n", j)) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* compute memory and flop estimates */
    /* ---------------------------------------------------------------------- */

    nchains = 0 ;			/* number of chains */
    Chain_start [0] = 0 ;		/* front 0 starts a new chain */
    maxrows = 1 ;
    maxcols = 1 ;
    maxfrsize = 1 ;
    dmaxfrsize = 1 ;
    chain_npiv = 0 ;			/* number of pivots in current chain */

    /* ---------------------------------------------------------------------- */
    /* simulate UMF_kernel_init */
    /* ---------------------------------------------------------------------- */

    /* Numeric->Memory usage estimate, in Units */
    head_usage = 1 ;		/* head marker (see UMF_mem_init_memoryspace) */
    init_tail_usage = 2 ;	/* tail marker (see UMF_mem_init_memoryspace) */
    dusage = 3 ;		/* head and tail markers */
    dhead_usage = 1 ;

    /* elements and tuples at tail*/
    init_tail_usage += UMF_kernel_init_usage (Ap, Row_degree, n_row, n_col,
	&dusage) ;

    /* ] done using Row_degree */

    ASSERT (UMF_is_permutation (Rperm_init, SWork->Ci, n_row, n_row)) ;

    tail_usage = init_tail_usage ;
    DEBUG2 (("tail_usage: "ID" (initial)\n", tail_usage)) ;
    max_usage = head_usage + init_tail_usage ;
    dmax_usage = dusage ;

    Symbolic->num_mem_init_usage = max_usage ;

    too_large = INT_OVERFLOW (dusage * sizeof (Unit)) ;
    if (too_large)
    {
	/* Initial memory usage, for input matrix only, */
	/* has encountered integer overflow. This is an error */
	/* condition, but keep going to compute other statistics. */
	Info [UMFPACK_VARIABLE_INIT_ESTIMATE] = dusage ;	/* too large */
    }
    else
    {
	Info [UMFPACK_VARIABLE_INIT_ESTIMATE] = (double) max_usage ;
    }

    /* ---------------------------------------------------------------------- */
    /* simulate UMF_kernel */
    /* ---------------------------------------------------------------------- */

    /* Use SWork->Ci as temporary workspace for link lists */
    Link = SWork->Ci ;
    for (i = 0 ; i < nfr ; i++)
    {
	Link [i] = EMPTY ;
    }

    dlnz = MIN (n_row, n_col) ;	/* upper limit of nz in L (incl diag) */
    dunz = dlnz ;		/* upper limit of nz in U (incl diag) */
    flops = 0 ;			/* flop count upper bound */

    DEBUG1 (("Umfpack symbolic, fronts:  nfr = "ID"\n", nfr)) ;

    for (i = 0 ; i < nfr ; i++)
    {

	fpivcol  = Front_npivcol [i] ; /* # candidate pivot columns */
	fallrows = SWork->Front_nrows [i] ;   /* all rows (not just Schur comp*/
	fallcols = SWork->Front_ncols [i] ;   /* all cols (not just Schur comp*/
	parent = Front_parent [i] ; /* parent in column etree */

	DEBUG1 (("\nFront "ID" fpivcol "ID" fallrows "ID" fallcols "ID"\n",
	    i, fpivcol, fallrows, fallcols)) ;

	/* determine the max size of the contribution block */
	fpiv = MIN (fpivcol, fallrows) ;	/* # pivot rows and cols */
	frows = fallrows - fpiv ;		/* max # rows in Schur comp. */
	fcols = fallcols - fpiv ;		/* max # cols in Schur comp. */

	DEBUG1 (("    "ID" : fpiv "ID" frows "ID" fcols "ID" parent "ID"\n",
		i, fpiv, frows, fcols, parent)) ;

	/* UMF_analyze can generate 0-by-c sized frontal matrices */
	ASSERT (fpiv >= 0 && fcols >= 0 && frows >= 0) ;

	/* assemble all children of front i in column etree */
	for (child = Link [i] ; child != EMPTY ; child = Link [child])
	{
	    ASSERT (child >= 0 && child < i) ;
	    ASSERT (Front_parent [child] == i) ;
	    /* free the child element */
	    cp = MIN (Front_npivcol [child], SWork->Front_nrows [child]);
	    cr = SWork->Front_nrows [child] - cp ;
	    cc = SWork->Front_ncols [child] - cp ;
	    ASSERT (cp >= 0 && cr >= 0 && cc >= 0) ;
	    tail_usage -= GET_ELEMENT_SIZE (cr, cc) + 1 ;
	    dusage -= DGET_ELEMENT_SIZE (cr, cc) + 1 ;
	    /* remove it from tuple lists */
	    tail_usage -= (cr + cc) * UNITS (Tuple, 1) ;
	    dusage -= ((double) cr + (double) cc) * UNITS (Tuple, 1) ;
	    DEBUG2 (("tail_usage: "ID" (assembled "ID" of size "ID")\n",
		tail_usage, child,
		GET_ELEMENT_SIZE (cr, cc) + 1 + (cr + cc) * UNITS (Tuple, 1))) ;
	}

	/* the flop count excludes the BLAS2 calls during pivot search */
	/* the assembly required to search a column, and assembly between */
	/* frontal matrices.  Those are "mushy" flops. */
	/* The flop count computed here is "canonical". */

	/* factorize the frontal matrix */
	f = (double) fpiv ;			/* # of pivots */
	r = (double) frows ;			/* # rows in Schur complement */
	c = (double) fcols ;			/* # cols in Schur complement */
	flops += DIV_FLOPS * (f*r + (f-1)*f/2)	/* scale pivot columns */
	    /* f outer products: */
	    + MULTSUB_FLOPS * (f*r*c + (r+c)*(f-1)*f/2 + (f-1)*f*(2*f-1)/6) ;

	/* count nonzeros in L and U */
	lf = (fpiv*fpiv-fpiv)/2 + fpiv*frows ;	/* nz in L below diagonal */
	uf = (fpiv*fpiv-fpiv)/2 + fpiv*fcols ;	/* nz in U above diagonal */

	/* store the f columns of L and f rows of U */
	head_usage +=
	    UNITS (Entry, lf + uf)	 /* numerical values (excl diagonal) */
	    + UNITS (Int, frows + fcols + fpiv) ; /* indices (compressed) */

	/* count nonzeros and memory usage in double precision */
	dlf = (f*f-f)/2 + f*r ;	/* nz in L below diagonal */
	duf = (f*f-f)/2 + f*c ;	/* nz in U above diagonal */
	dlnz += dlf ;
	dunz += duf ;
	dh =
	    DUNITS (Entry, dlf + duf)    /* numerical values (excl diagonal) */
	    + DUNITS (Int, r + c + f) ; /* indices (compressed) */
	dusage += dh ;
	dhead_usage += dh ;

	if (parent != EMPTY)
	{
	    /* create new element */
	    tail_usage += GET_ELEMENT_SIZE (frows, fcols) + 1 ;
	    dusage += DGET_ELEMENT_SIZE (frows, fcols) + 1 ;

	    /* place new element in tuple lists */
	    tail_usage += (frows + fcols) * UNITS (Tuple, 1) ;
	    dusage += (r + c) * UNITS (Tuple, 1) ;
	    DEBUG2 (("tail_usage: "ID" (create    "ID" of size "ID")\n",
		tail_usage, i,
		GET_ELEMENT_SIZE (frows, fcols) + 1
		+ (frows + fcols) * UNITS (Tuple, 1))) ;

	    /* place in link list of parent */
	    Link [i] = Link [parent] ;
	    Link [parent] = i ;
	}

	/* keep track of peak Numeric->Memory usage */
	max_usage = MAX (max_usage, head_usage + tail_usage) ;

	/* max_usage may encounter integer overflow, so dmax_usage also kept. */
	/* account for possible roundoff errors in dusage. */
	/* Ignore NaN case. */
	dusage *= (1.0 + MAX_EPSILON) ;
	dusage = MAX (dusage,
	    (double) (head_usage + tail_usage) * (1.0 + MAX_EPSILON)) ;
	dmax_usage = MAX (dmax_usage, dusage) ;
	dhead_usage *= (1.0 + MAX_EPSILON) ;
	dhead_usage = MAX (dhead_usage,
	    (double) head_usage * (1.0 + MAX_EPSILON)) ;

	/* at most nb or chain_npiv pending pivots from old fronts */
	maxoldpiv = MIN (nb, chain_npiv) ;
	maxrows = MAX (maxrows, maxoldpiv + fallrows) ;
	maxcols = MAX (maxcols, maxoldpiv + fallcols) ;

	DEBUG1 (("    ID: chain_npiv "ID" fpiv "ID" fallrows "ID,
		i, chain_npiv, fpiv, fallrows)) ;
	DEBUG1 ((" fallcols "ID"   maxoldpiv "ID" maxrows "ID" maxcols "ID"\n",
		fallcols, maxoldpiv, maxrows, maxcols)) ;

	chain_npiv += fpiv ;
	if (parent != i+1)
	{
	    /* this is the end of a chain */
	    Chain_maxrows [nchains] = maxrows ;
	    Chain_maxcols [nchains] = maxcols ;
	    /* Ignore NaN case. */
	    dmaxfrsize = MAX (dmaxfrsize, (double) maxrows * (double) maxcols) ;
	    maxfrsize = MAX (maxfrsize, maxrows * maxcols) ;
	    nchains++ ;
	    Chain_start [nchains] = i+1 ;
	    maxrows = 1 ;
	    maxcols = 1 ;
	    chain_npiv = 0 ;
	}
    }

    dhead_usage = ceil (dhead_usage) ;
    dmax_usage = ceil (dmax_usage) ;

    tail_usage -= init_tail_usage ;

    /* all tuples and elements are now deallocated */
    DEBUG0 (("final tail_usage: "ID"\n", tail_usage)) ;
    ASSERT (tail_usage == 0) ;

    DEBUG1 (("dmaxfrsize %30.20g Int_MAX %30d\n", dmaxfrsize, Int_MAX)) ;

    /* check if the frontal matrix is too big */
    too_large = too_large || INT_OVERFLOW (dmaxfrsize * sizeof (Entry)) ;

    /* ---------------------------------------------------------------------- */
    /* find the biggest frontal matrix, for all chains */
    /* ---------------------------------------------------------------------- */

    maxnrows = 1 ;
    maxncols = 1 ;
    for (i = 0 ; i < nchains ; i++)
    {
	maxnrows = MAX (maxnrows, Chain_maxrows [i]) ;
	maxncols = MAX (maxncols, Chain_maxcols [i]) ;
    }

    /* information to keep for numeric factorization */
    Symbolic->maxfrsize = maxfrsize ;
    Symbolic->maxnrows = maxnrows ;
    Symbolic->maxncols = maxncols ;
    Symbolic->num_mem_usage_est = dmax_usage ;
    Symbolic->num_mem_size_est = dhead_usage ;

    /* ---------------------------------------------------------------------- */
    /* estimate total memory usage in UMFPACK_numeric */
    /* ---------------------------------------------------------------------- */

    UMF_set_stats (
	Info,
	Symbolic,
	dmax_usage,		/* estimated peak size of Numeric->Memory */
	dhead_usage,		/* estimated final size of Numeric->Memory */
	flops,			/* estimated "true flops" */
	dlnz,			/* estimated nz in L */
	dunz,			/* estimated nz in U */
	dmaxfrsize,		/* estimated largest front size */
	(double) n_col,		/* worst case Numeric->Upattern size */
	(double) n_inner,	/* max possible pivots to be found */
	ESTIMATE) ;

    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    for (i = 0 ; i < nchains ; i++)
    {
	DEBUG2 (("Chain "ID" start "ID" end "ID" maxrows "ID" maxcols "ID"\n",
		i, Chain_start [i], Chain_start [i+1] - 1,
		Chain_maxrows [i], Chain_maxcols [i])) ;
	UMF_dump_chain (Chain_start [i],
	    SWork->Front_parent,
	    SWork->Front_npivcol,
	    SWork->Front_nrows,
	    SWork->Front_ncols,
	    nfr) ;
    }
    fpivcol = 0 ;
    for (i = 0 ; i < nfr ; i++)
    {
	fpivcol = MAX (fpivcol, Front_npivcol [i]) ;
    }
    DEBUG0 (("Max pivot cols in any front: "ID"\n", fpivcol)) ;
    DEBUG1 (("Largest front: maxnrows "ID" maxncols "ID" maxfrsize "ID"\n",
	maxnrows, maxncols, maxfrsize)) ;
    DEBUG1 (("(savings "ID") nchains "ID"\n",
	(maxnrows*maxncols) - maxfrsize, nchains)) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* is the problem too large? */
    /* ---------------------------------------------------------------------- */

    if (too_large)
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_problem_too_large ;
	error (&Symbolic, SWork) ;
	return (UMFPACK_ERROR_problem_too_large) ;
    }

    /* ---------------------------------------------------------------------- */
    /* UMFPACK_symbolic was successful, return the object handle */
    /* ---------------------------------------------------------------------- */

    Symbolic->valid = SYMBOLIC_VALID ;
    *SymbolicHandle = (void *) Symbolic ;

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    free_work (SWork) ;
    /* Symbolic contains 10 objects */
    DEBUG0 (("(3)Symbolic UMF_malloc_count - init_count = "ID"\n",
	UMF_malloc_count - init_count)) ;
    ASSERT (UMF_malloc_count == init_count + 10) ;

    /* ---------------------------------------------------------------------- */
    /* get the time used by UMFPACK_*symbolic */
    /* ---------------------------------------------------------------------- */

    tend = umfpack_timer ( ) ;
    Info [UMFPACK_SYMBOLIC_TIME] = MAX (0, tend - tstart) ;

    return (UMFPACK_OK) ;
}


/* ========================================================================== */
/* === free_work ============================================================ */
/* ========================================================================== */

PRIVATE void free_work
(
    SWorkType *SWork
)
{
    ASSERT (SWork) ;

    SWork->Ci = (Int *) UMF_free ((void *) SWork->Ci) ;
    SWork->Front_npivcol = (Int *) UMF_free ((void *) SWork->Front_npivcol) ;
    SWork->Front_nrows = (Int *) UMF_free ((void *) SWork->Front_nrows) ;
    SWork->Front_ncols = (Int *) UMF_free ((void *) SWork->Front_ncols) ;
    SWork->Front_parent = (Int *) UMF_free ((void *) SWork->Front_parent) ;
    SWork->Front_cols = (Int *) UMF_free ((void *) SWork->Front_cols) ;

}


/* ========================================================================== */
/* === error ================================================================ */
/* ========================================================================== */

/* Error return from UMFPACK_symbolic.  Free all allocated memory. */

PRIVATE void error
(
    SymbolicType **Symbolic,
    SWorkType *SWork
)
{

    free_work (SWork) ;
    UMFPACK_free_symbolic ((void **) Symbolic) ;
    ASSERT (UMF_malloc_count == init_count) ;
}

