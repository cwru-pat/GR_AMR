/* ========================================================================== */
/* === UMFPACK_numeric ====================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Factorizes A into its LU factors, given a symbolic
    pre-analysis computed by UMFPACK_symbolic.  See umfpack_numeric.h for a
    description.

    Dynamic memory usage: UMFPACK_numeric makes the heaviest usage of memory
    space of all UMFPACK routines.  Here is an outline of how it uses memory:

	1) calls UMF_malloc 14 times, to obtain temporary workspace of size f
	   Entry's and 2*(n_row+1) + 2*(n_col+1) + (n_col+n_inner+1) +
	   (nn+1) + 3*(c+1) + 2*(r+1) + max(r,c) + (nfr+1) integers
	   where f is the size of the working array used for the frontal
	   matrices, r is the maximum number of rows in the working array,
	   c is the maximum number of columns in the working array (f is
	   less than or equal to r*c, and typically f=r*c), n_inner is
	   min (n_row,n_col), nn is max (n_row,n_col), and nfr is the number
	   of frontal matrices.  For a square matrix, this is f Entry's and
	   about 7n + 3c + 2r + max(r,c) + nfr integers.

	2) calls UMF_malloc 10 times, for a total space of sizeof (NumericType)
	   bytes, 4*(n_row+1) + 4*(n_row+1) integers, and (n_inner+1) Entry's.
	   sizeof (NumericType) is a small constant.  This space is part of the
	   permanent Numeric object (which holds the LU factors) and is not
	   freed on return via UMF_free unless an error occurs.

	3) calls UMF_malloc once, to allocate the variable-sized part of the
	   Numeric object, whose size in Units is the larger of:
	   (Control [UMFPACK_ALLOC_INIT]) *  (the approximate upper bound
	   computed by UMFPACK_symbolic), and the minimum required to start the
	   numerical factorization.  The default value of
	   Control [UMFPACK_ALLOC_INIT] is 1.0.  This request is reduced if it
	   fails.  This object is not freed on return via UMF_free unless
	   an error occurs.

	4) During numerical factorization (inside UMF_kernel), the variable-size
	   block of memory is increased in size via a call to UMF_realloc if
	   it is found to be too small.  This is rare with the default control
	   setting of 1.0 (I've never observed it to happen, but it
	   theoretically could happen).  During factorization, this block holds
	   the pattern and values of L and U at the top end, and the elements
	   (contibution blocks) at the bottom end.

	   For a square nonsingular matrix, the peak memory usage at this point
	   is roughly:

		r*c + n Entry's
		15n + 3c + 2r + max(r,c) + nfr integers
		plus peak size of the variable-sized part of the Numeric object
		plus the size of the Symbolic object (2n to 9n integers)

	   where r is the maximum number of rows in the working array used
	   to hold the current frontal matrix, and c is the maximum number of
	   columns in the working array.

	   The peak value of the variable-sized object is estimated in
	   UMFPACK_*symbolic (Info [UMFPACK_VARIABLE_PEAK_ESTIMATE]).
	   The size of the Symbolic object is in Info [UMFPACK_SYMBOLIC_SIZE],
	   and is between 2*n and 9*n integers.

	5) After numerical factorization all of the objects allocated in step
	   (1) are freed via UMF_free, except that one object of size n_col+1
	   is kept if there are nonzeros in the last pivot row.

	6) The variable-sized block is reduced to hold just L and U, via a call
	   to UMF_realloc, since the frontal matrices are no longer needed.

	7) This leaves a total of 11 or 12 objects allocated by UMF_malloc that
	   form the LU factorization.  These remain if UMFPACK_numeric was
	   successful.  Otherwise, they are all freed via UMF_free.
	   The final size of the Numeric object for a square nonsingular matrix
	   is roughly:

		n Entry's
		8n integers
		plus final size of the variable-sized part of the Numeric object

    Dynamic memory usage of UMFPACK_free_numeric:

	1) It frees, via UMF_free, all 11 or 12 objects allocated for the
	   Numeric object, in a prior call to UMFPACK_numeric.

*/

/* ========================================================================== */

#include "umf_internal.h"
#include "umf_valid_symbolic.h"
#include "umf_set_stats.h"
#include "umf_kernel.h"
#include "umf_malloc.h"
#include "umf_free.h"
#include "umf_realloc.h"

#ifndef NDEBUG
PRIVATE Int init_count ;
#endif

PRIVATE Int work_alloc
(
    WorkType *Work
) ;

PRIVATE void free_work
(
    WorkType *Work
) ;

PRIVATE Int numeric_alloc
(
    NumericType **NumericHandle,
    SymbolicType *Symbolic,
    double alloc_init
) ;

PRIVATE void error
(
    NumericType **Numeric,
    WorkType *Work
) ;


/* ========================================================================== */

GLOBAL Int UMFPACK_numeric
(
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    void *SymbolicHandle,
    void **NumericHandle,
    const double Control [UMFPACK_CONTROL],
    double User_Info [UMFPACK_INFO]
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    NumericType *Numeric ;
    SymbolicType *Symbolic ;
    WorkType WorkSpace, *Work ;
    Int n_row, n_col, n_inner, newsize, i, status, *inew, npiv, ulen ;
    Unit *mnew ;
    double Info2 [UMFPACK_INFO], *Info, alloc_init, relax, relpt, tstart, tend,
	relax2, relax3 ;

    /* ---------------------------------------------------------------------- */
    /* get the amount of time used by the process so far */
    /* ---------------------------------------------------------------------- */

    tstart = umfpack_timer ( ) ;

    /* ---------------------------------------------------------------------- */
    /* initialize and check inputs */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    init_count = UMF_malloc_count ;
#endif

    if (Control)
    {
	/* use the Control array passed to us by the caller; look for NaN's */

	if (SCALAR_IS_NAN (Control [UMFPACK_PIVOT_TOLERANCE]))
	{
	    relpt = UMFPACK_DEFAULT_PIVOT_TOLERANCE ;
	}
	else
	{
	    relpt = Control [UMFPACK_PIVOT_TOLERANCE] ;
	}

	if (SCALAR_IS_NAN (Control [UMFPACK_RELAXED_AMALGAMATION]))
	{
	    relax = UMFPACK_DEFAULT_RELAXED_AMALGAMATION ;
	}
	else
	{
	    relax = Control [UMFPACK_RELAXED_AMALGAMATION] ;
	}

	if (SCALAR_IS_NAN (Control [UMFPACK_RELAXED2_AMALGAMATION]))
	{
	    relax2 = UMFPACK_DEFAULT_RELAXED2_AMALGAMATION ;
	}
	else
	{
	    relax2 = Control [UMFPACK_RELAXED2_AMALGAMATION] ;
	}

	if (SCALAR_IS_NAN (Control [UMFPACK_RELAXED3_AMALGAMATION]))
	{
	    relax3 = UMFPACK_DEFAULT_RELAXED3_AMALGAMATION ;
	}
	else
	{
	    relax3 = Control [UMFPACK_RELAXED3_AMALGAMATION] ;
	}

	if (SCALAR_IS_NAN (Control [UMFPACK_ALLOC_INIT]))
	{
	    alloc_init = UMFPACK_DEFAULT_ALLOC_INIT ;
	}
	else
	{
	    alloc_init = Control [UMFPACK_ALLOC_INIT] ;
	}

    }
    else
    {
	/* no Control passed - use defaults instead */
	relpt = UMFPACK_DEFAULT_PIVOT_TOLERANCE ;
	relax = UMFPACK_DEFAULT_RELAXED_AMALGAMATION ;
	relax2 = UMFPACK_DEFAULT_RELAXED2_AMALGAMATION ;
	relax3 = UMFPACK_DEFAULT_RELAXED3_AMALGAMATION ;
	alloc_init = UMFPACK_DEFAULT_ALLOC_INIT ;
    }

    relpt = MAX (0.0, MIN (relpt, 1.0)) ;
    relax = MAX (0.0, relax) ;
    relax2 = MAX (0.0, relax2) ;
    relax3 = MAX (0.0, relax3) ;
    alloc_init = MAX (0.0, alloc_init) ;

    if (User_Info)
    {
	/* return Info in user's array */
	Info = User_Info ;
	for (i = UMFPACK_NUMERIC_SIZE ; i <= UMFPACK_NUMERIC_TIME ; i++)
	{
	    Info [i] = EMPTY ;
	}
    }
    else
    {
	/* no Info array passed - use local one instead */
	Info = Info2 ;
	for (i = 0 ; i < UMFPACK_INFO ; i++)
	{
	    Info [i] = EMPTY ;
	}
    }

    Symbolic = (SymbolicType *) SymbolicHandle ;
    Numeric = (NumericType *) NULL ;
    if (!UMF_valid_symbolic (Symbolic))
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_invalid_Symbolic_object ;
	return (UMFPACK_ERROR_invalid_Symbolic_object) ;
    }

    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;
    n_inner = MIN (n_row, n_col) ;

    Info [UMFPACK_STATUS] = UMFPACK_OK ;
    Info [UMFPACK_NROW] = n_row ;
    Info [UMFPACK_NCOL] = n_col ;
    Info [UMFPACK_SIZE_OF_UNIT] = (double) (sizeof (Unit)) ;
    Info [UMFPACK_NUMERIC_DEFRAG] = 0 ;
    Info [UMFPACK_NUMERIC_REALLOC] = 0 ;
    Info [UMFPACK_NUMERIC_COSTLY_REALLOC] = 0 ;

    if (!Ap || !Ai || !Ax || !NumericHandle
#ifdef COMPLEX
	|| !Az
#endif
    )
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_argument_missing ;
	return (UMFPACK_ERROR_argument_missing) ;
    }

    Info [UMFPACK_NZ] = Ap [n_col] ;
    *NumericHandle = (void *) NULL ;

    /* ---------------------------------------------------------------------- */
    /* allocate the Work object */
    /* ---------------------------------------------------------------------- */

    Work = &WorkSpace ;
    Work->n_row = n_row ;
    Work->n_col = n_col ;
    Work->nfr = Symbolic->nfr ;
    Work->maxfrsize = Symbolic->maxfrsize ;
    Work->maxnrows = Symbolic->maxnrows ;
    Work->maxncols = Symbolic->maxncols ;

    if (!work_alloc (Work))
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	error (&Numeric, Work) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }
    ASSERT (UMF_malloc_count == init_count + 14) ;

    /* ---------------------------------------------------------------------- */
    /* allocate Numeric object */
    /* ---------------------------------------------------------------------- */

    if (!numeric_alloc (&Numeric, Symbolic, alloc_init))
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	error (&Numeric, Work) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }
    DEBUG0 (("malloc: init_count "ID" UMF_malloc_count "ID"\n",
	init_count, UMF_malloc_count)) ;
    ASSERT (UMF_malloc_count == init_count + 14 + 11) ;

    /* set control parameters */
    Numeric->relpt = relpt ;
    Numeric->relax = relax ;
    Numeric->relax2 = relax2 ;
    Numeric->relax3 = relax3 ;
    Numeric->alloc_init = alloc_init ;

    DEBUG0 (("umf control: relpt %g relax [%g %g %g] init %g inc %g red %g\n",
	relpt, relax, relax2, relax3, alloc_init,
	UMF_REALLOC_INCREASE, UMF_REALLOC_REDUCTION)) ;

    /* ---------------------------------------------------------------------- */
    /* factorize */
    /* ---------------------------------------------------------------------- */

    status = UMF_kernel (Ap, Ai, Ax,
#ifdef COMPLEX
	Az,
#endif
	Numeric, Work, Symbolic) ;
    Info [UMFPACK_STATUS] = status ;
    Info [UMFPACK_VARIABLE_INIT] = Numeric->init_usage ;
    if (status < 0)
    {
	/* out of memory, or pattern has changed */
	error (&Numeric, Work) ;
	return (status) ;
    }
    DEBUG0 (("malloc: init_count "ID" UMF_malloc_count "ID"\n",
	init_count, UMF_malloc_count)) ;

    npiv = Numeric->npiv ;	/* = n_inner for nonsingular matrices */
    ulen = Numeric->ulen ;	/* = 0 for square nonsingular matrices */

    /* ---------------------------------------------------------------------- */
    /* free Work object */
    /* ---------------------------------------------------------------------- */

    DEBUG0 (("malloc: init_count "ID" UMF_malloc_count "ID"\n",
	init_count, UMF_malloc_count)) ;
    free_work (Work) ;
    DEBUG0 (("malloc: init_count "ID" UMF_malloc_count "ID"\n",
	init_count, UMF_malloc_count)) ;
    DEBUG0 (("Numeric->ulen: "ID"\n", ulen)) ;
    ASSERT (UMF_malloc_count == init_count + 11 + (ulen > 0)) ;

    /* ---------------------------------------------------------------------- */
    /* reduce Lpos, Lilen, Lip, Upos, Uilen and Uip to size npiv+1 */
    /* ---------------------------------------------------------------------- */

    /* This is only needed for rectangular or singular matrices. */

    if (npiv < n_row)
    {
	/* reduce Lpos, Uilen, and Uip from size n_row+1 to size npiv */
	inew = (Int *) UMF_realloc (Numeric->Lpos, npiv+1, sizeof (Int)) ;
	if (inew)
	{
	    Numeric->Lpos = inew ;
	}
	inew = (Int *) UMF_realloc (Numeric->Uilen, npiv+1, sizeof (Int)) ;
	if (inew)
	{
	    Numeric->Uilen = inew ;
	}
	inew = (Int *) UMF_realloc (Numeric->Uip, npiv+1, sizeof (Int)) ;
	if (inew)
	{
	    Numeric->Uip = inew ;
	}
    }

    if (npiv < n_col)
    {
	/* reduce Upos, Lilen, and Lip from size n_col+1 to size npiv */
	inew = (Int *) UMF_realloc (Numeric->Upos, npiv+1, sizeof (Int)) ;
	if (inew)
	{
	    Numeric->Upos = inew ;
	}
	inew = (Int *) UMF_realloc (Numeric->Lilen, npiv+1, sizeof (Int)) ;
	if (inew)
	{
	    Numeric->Lilen = inew ;
	}
	inew = (Int *) UMF_realloc (Numeric->Lip, npiv+1, sizeof (Int)) ;
	if (inew)
	{
	    Numeric->Lip = inew ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* reduce Numeric->Upattern from size n_col+1 to size ulen+1 */
    /* ---------------------------------------------------------------------- */

    /* This is only needed for singular matrices, or if n_row < n_col. */
    /* If ulen is zero, the object does not exist. */

    DEBUG4 (("ulen: "ID" Upattern "ID"\n", ulen, (Int) Numeric->Upattern)) ;
    ASSERT (IMPLIES (ulen == 0, Numeric->Upattern == (Int *) NULL)) ;
    if (ulen > 0 && ulen < n_col)
    {
	inew = (Int *) UMF_realloc (Numeric->Upattern, ulen+1, sizeof (Int)) ;
	if (inew)
	{
	    Numeric->Upattern = inew ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* reduce Numeric->Memory to hold just the LU factors at the head */
    /* ---------------------------------------------------------------------- */

    newsize = Numeric->ihead ;
    if (newsize < Numeric->size)
    {
	mnew = (Unit *) UMF_realloc (Numeric->Memory, newsize, sizeof (Unit)) ;
	if (mnew)
	{
	    /* realloc succeeded (how can it fail since the size is reduced?) */
	    Numeric->Memory = mnew ;
	    Numeric->size = newsize ;
	}
    }
    Numeric->ihead = Numeric->size ;
    Numeric->itail = Numeric->ihead ;
    Numeric->tail_usage = 0 ;
    Numeric->ibig = EMPTY ;
    /* UMF_mem_alloc_tail_block can no longer be called (no tail marker) */

    /* ---------------------------------------------------------------------- */
    /* report the results and return the Numeric object */
    /* ---------------------------------------------------------------------- */

    UMF_set_stats (
	Info,
	Symbolic,
	(double) Numeric->max_usage,	/* actual peak Numeric->Memory */
	(double) Numeric->size,		/* actual final Numeric->Memory */
	Numeric->flops,			/* actual "true flops" */
	(double) Numeric->lnz + n_inner,		/* actual nz in L */
	(double) Numeric->unz + Numeric->nnzpiv,	/* actual nz in U */
	(double) Numeric->maxfrsize,	/* actual largest front size */
	(double) ulen,			/* actual Numeric->Upattern size */
	(double) npiv,			/* actual # pivots found */
	ACTUAL) ;

    Info [UMFPACK_NUMERIC_DEFRAG] = Numeric->ngarbage ;
    Info [UMFPACK_NUMERIC_REALLOC] = Numeric->nrealloc ;
    Info [UMFPACK_NUMERIC_COSTLY_REALLOC] = Numeric->ncostly ;
    Info [UMFPACK_COMPRESSED_PATTERN] = Numeric->isize ;
    Info [UMFPACK_LU_ENTRIES] = Numeric->nLentries + Numeric->nUentries +
    	Numeric->npiv ;
    Info [UMFPACK_UDIAG_NZ] = Numeric->nnzpiv ;

    /* estimate of the recipricol of the condition number. */
    if (SCALAR_IS_ZERO (Numeric->min_udiag)
     || SCALAR_IS_ZERO (Numeric->max_udiag))
    {
	/* rcond is zero if there is any zero on the diagonal, */
	/* even if NaN's are also present. */
	Numeric->rcond = 0.0 ;
    }
    else
    {
	/* estimate of the recipricol of the condition number. */
	/* This is NaN if diagonal is zero-free, but has one or more NaN's. */
	Numeric->rcond = Numeric->min_udiag / Numeric->max_udiag ;
    }
    Info [UMFPACK_RCOND] = Numeric->rcond ;

    if (Numeric->nnzpiv < n_inner
    || SCALAR_IS_ZERO (Numeric->rcond) || SCALAR_IS_NAN (Numeric->rcond))
    {
	/* there are zeros and/or NaN's on the diagonal of U */
	DEBUG0 (("Warning, matrix is singular in umfpack_numeric\n")) ;
	DEBUG0 (("nnzpiv "ID" n_inner "ID" rcond %g\n", Numeric->nnzpiv,
	    n_inner, Numeric->rcond)) ;
	status = UMFPACK_WARNING_singular_matrix ;
	Info [UMFPACK_STATUS] = status ;
    }

    Numeric->valid = NUMERIC_VALID ;
    *NumericHandle = (void *) Numeric ;

    /* Numeric has 11 or 12 objects */
    ASSERT (UMF_malloc_count == init_count + 11 + (ulen > 0)) ;

    /* ---------------------------------------------------------------------- */
    /* get the time used by UMFPACK_numeric */
    /* ---------------------------------------------------------------------- */

    tend = umfpack_timer ( ) ;
    Info [UMFPACK_NUMERIC_TIME] = MAX (0, tend - tstart) ;

    /* return UMFPACK_OK or UMFPACK_WARNING_singular_matrix */
    return (status) ;

}


/* ========================================================================== */
/* === numeric_alloc ======================================================== */
/* ========================================================================== */

/* Allocate the Numeric object */

PRIVATE Int numeric_alloc
(
    NumericType **NumericHandle,
    SymbolicType *Symbolic,
    double alloc_init
)
{
    Int n_row, n_col, n_inner, min_usage, trying ;
    NumericType *Numeric ;
    double nsize, bsize ;

    ASSERT (Symbolic) ;

    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;
    n_inner = MIN (n_row, n_col) ;
    *NumericHandle = (NumericType *) NULL ;

    /* 1 allocation:  accounted for in UMF_set_stats (num_fixed_size) */
    Numeric = (NumericType *) UMF_malloc (1, sizeof (NumericType)) ;

    if (!Numeric)
    {
	return (FALSE) ;	/* out of memory */
    }
    Numeric->valid = 0 ;
    *NumericHandle = Numeric ;

    /* 9 allocations:  accounted for in UMF_set_stats (num_fixed_size) */
    Numeric->D = (Entry *) UMF_malloc (n_inner+1, sizeof (Entry)) ;
    Numeric->Rperm = (Int *) UMF_malloc (n_row+1, sizeof (Int)) ;
    Numeric->Cperm = (Int *) UMF_malloc (n_col+1, sizeof (Int)) ;
    Numeric->Lpos = (Int *) UMF_malloc (n_row+1, sizeof (Int)) ;
    Numeric->Lilen = (Int *) UMF_malloc (n_col+1, sizeof (Int)) ;
    Numeric->Lip = (Int *) UMF_malloc (n_col+1, sizeof (Int)) ;
    Numeric->Upos = (Int *) UMF_malloc (n_col+1, sizeof (Int)) ;
    Numeric->Uilen = (Int *) UMF_malloc (n_row+1, sizeof (Int)) ;
    Numeric->Uip = (Int *) UMF_malloc (n_row+1, sizeof (Int)) ;

    Numeric->Memory = (Unit *) NULL ;
    Numeric->Upattern = (Int *) NULL ;	/* used for singular matrices only */

    if (!Numeric->D || !Numeric->Rperm || !Numeric->Cperm || !Numeric->Upos ||
	!Numeric->Lpos || !Numeric->Lilen || !Numeric->Uilen || !Numeric->Lip ||
	!Numeric->Uip)
    {
	return (FALSE) ;	/* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* allocate initial Numeric->Memory for LU factors and elements */
    /* ---------------------------------------------------------------------- */

    nsize = (alloc_init * Symbolic->num_mem_usage_est) * (1.0 + MAX_EPSILON) ;
    nsize = ceil (nsize) ;
    min_usage = Symbolic->num_mem_init_usage ;

    /* Numeric->Memory must be large enough for UMF_kernel_init */
    /* double relop, but ignore NaN case. */
    nsize = MAX (ceil ((double) min_usage * (1.0 + MAX_EPSILON)), nsize) ;

    /* Numeric->Memory cannot be larger in size than Int_MAX / sizeof(Unit) */
    /* For ILP32 mode:  2GB (nsize cannot be bigger than 256 Mwords) */

    bsize = ((double) Int_MAX) / sizeof (Unit) - 1 ;
    DEBUG0 (("bsize %g\n", bsize)) ;

    nsize = MIN (nsize, bsize) ;    /* double relop, but ignore NaN case. */

    Numeric->size = (Int) nsize ;

    DEBUG0 (("Num init %g usage_est %g numsize "ID" minusage "ID"\n",
	alloc_init, Symbolic->num_mem_usage_est, Numeric->size, min_usage)) ;

    /* allocates 1 object: */
    /* keep trying until successful, or memory request is too small */
    trying = TRUE ;
    while (trying)
    {
	Numeric->Memory = (Unit *) UMF_malloc (Numeric->size, sizeof (Unit)) ;
	if (Numeric->Memory)
	{
	    DEBUG0 (("Successful Numeric->size: "ID"\n", Numeric->size)) ;
	    return (TRUE) ;
	}
	/* too much, reduce the request (but not below the minimum) */
	/* and try again */
	trying = Numeric->size > min_usage ;
	Numeric->size = (Int)
	    (UMF_REALLOC_REDUCTION * ((double) Numeric->size)) ;
	Numeric->size = MAX (min_usage, Numeric->size) ;
    }

    return (FALSE) ;	/* we failed to allocate Numeric->Memory */
}


/* ========================================================================== */
/* === work_alloc =========================================================== */
/* ========================================================================== */

/* Allocate the Work object.  Return TRUE if successful. */

PRIVATE Int work_alloc
(
    WorkType *Work
)
{
    Int n_row, n_col, nn, n_inner, maxfrsize, maxnrows, maxncols, nfr ;

    n_row = Work->n_row ;
    n_col = Work->n_col ;
    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;
    nfr = Work->nfr ;
    maxfrsize = Work->maxfrsize ;
    maxnrows = Work->maxnrows ;
    maxncols = Work->maxncols ;

    /* largest front will be maxnrows-by-maxncols, at most */
    DEBUG1 (("Allocating frontal matrix, size "ID" ("ID"-by-"ID") saved: "ID
	"\n", maxfrsize, maxnrows, maxncols,
	(maxnrows * maxncols) - maxfrsize)) ;

    /* 13 allocations, freed in free_work: */
    /* accounted for in UMF_set_stats (work_usage) */
    Work->Fx = (Entry *) UMF_malloc (maxfrsize, sizeof (Entry)) ;
    Work->Frpos = (Int *) UMF_malloc (n_row + 1, sizeof (Int)) ;
    Work->Fcpos = (Int *) UMF_malloc (n_col + 1, sizeof (Int)) ;
    Work->Lpattern = (Int *) UMF_malloc (n_row + 1, sizeof (Int)) ;
    Work->Wp = (Int *) UMF_malloc (nn + 1, sizeof (Int)) ;
    Work->Frows = (Int *) UMF_malloc (maxnrows + 1, sizeof (Int)) ;
    Work->Fcols = (Int *) UMF_malloc (maxncols + 1, sizeof (Int)) ;
    Work->Wio = (Int *) UMF_malloc (maxncols + 1, sizeof (Int)) ;
    Work->Woi = (Int *) UMF_malloc (maxncols + 1, sizeof (Int)) ;
    Work->Woo = (Int *) UMF_malloc (MAX (maxnrows, maxncols) + 1, sizeof (Int));
    Work->Wm = (Int *) UMF_malloc (maxnrows + 1, sizeof (Int)) ;
    Work->E = (Int *) UMF_malloc (n_col + n_inner + 1, sizeof (Int)) ;
    Work->Front_new1strow = (Int *) UMF_malloc (nfr + 1, sizeof (Int)) ;

    /* 1 allocation, may become part of Numeric (if singular or rectangular): */
    Work->Upattern = (Int *) UMF_malloc (n_col + 1, sizeof (Int)) ;

    return (Work->Fx && Work->Frpos
	&& Work->Fcpos && Work->Lpattern && Work->Upattern
	&& Work->E && Work->Frows && Work->Fcols && Work->Wio
	&& Work->Woi && Work->Woo && Work->Wm && Work->Wp
	&& Work->Front_new1strow) ;
}


/* ========================================================================== */
/* === free_work ============================================================ */
/* ========================================================================== */

PRIVATE void free_work
(
    WorkType *Work
)
{
    if (Work)
    {

	/* free 13 or 14 objects (Upattern may already be gone) */
	Work->Fx = (Entry *) UMF_free ((void *) Work->Fx) ;
	Work->Frpos = (Int *) UMF_free ((void *) Work->Frpos) ;
	Work->Fcpos = (Int *) UMF_free ((void *) Work->Fcpos) ;
	Work->Lpattern = (Int *) UMF_free ((void *) Work->Lpattern) ;
	Work->Upattern = (Int *) UMF_free ((void *) Work->Upattern) ;
	Work->Wp = (Int *) UMF_free ((void *) Work->Wp) ;
	Work->Frows = (Int *) UMF_free ((void *) Work->Frows) ;
	Work->Fcols = (Int *) UMF_free ((void *) Work->Fcols) ;
	Work->Wio = (Int *) UMF_free ((void *) Work->Wio) ;
	Work->Woi = (Int *) UMF_free ((void *) Work->Woi) ;
	Work->Woo = (Int *) UMF_free ((void *) Work->Woo) ;
	Work->Wm = (Int *) UMF_free ((void *) Work->Wm) ;
	Work->E = (Int *) UMF_free ((void *) Work->E) ;
	Work->Front_new1strow =
	    (Int *) UMF_free ((void *) Work->Front_new1strow) ;

    }
}


/* ========================================================================== */
/* === error ================================================================ */
/* ========================================================================== */

/* Error return from UMFPACK_numeric.  Free all allocated memory. */

PRIVATE void error
(
    NumericType **Numeric,
    WorkType *Work
)
{
    free_work (Work) ;
    UMFPACK_free_numeric ((void **) Numeric) ;
    ASSERT (UMF_malloc_count == init_count) ;
}

