/* ========================================================================== */
/* === UMF_kernel_init ====================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Initialize the kernel, build tuple lists.  Assumes elements are packed.
    Returns TRUE if successful, FALSE if out of memory or if the pattern has
    changed since UMFPACK_*symbolic.  UMFPACK_numeric allocates at least enough
    space for UMF_kernel_init to succeed; otherwise it does not call
    UMF_kernel_init.  So an out-of-memory condition means that the pattern must
    have gotten larger.
*/

#include "umf_internal.h"
#include "umf_tuple_lengths.h"
#include "umf_build_tuples.h"
#include "umf_mem_init_memoryspace.h"
#include "umf_mem_alloc_element.h"

GLOBAL Int UMF_kernel_init
(
    const Int Ap [ ],		/* user's input matrix (not modified) */
    const Int Ai [ ],
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int row, k, oldcol, size, e, p1, p2, p, *ip, nz, *Rows, *Cols, *E, i, *Upos,
	*Lpos, n_row, n_col, *Wp, *Cperm_init, *Frpos, *Fcpos, *Row_degree, nn,
	*Row_tlen, *Col_degree, *Col_tlen, deg, oldrow, newrow, ilast,
	*Rperm_init, col, n_inner ;
    double unused = 0 ;
    Entry *xp, *C ;
    Element *ep ;

#ifndef NDEBUG
    double gprob_save = UMF_gprob ;
    UMF_gprob = 0 ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    DEBUG0 (("KERNEL INIT\n")) ;

    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;
    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;
    ASSERT (n_row > 0 && n_col > 0) ;
    Cperm_init = Symbolic->Cperm_init ;
    Rperm_init = Symbolic->Rperm_init ;
    nz = Ap [n_col] ;
    if (nz < 0 || Ap [0] != 0 || nz != Symbolic->nz)
    {
	return (FALSE) ;	/* pattern changed */
    }

    /* ---------------------------------------------------------------------- */
    /* initialize the Numeric->Memory space for LU, elements, and tuples */
    /* ---------------------------------------------------------------------- */

    UMF_mem_init_memoryspace (Numeric) ;

    /* ---------------------------------------------------------------------- */
    /* initialize the Work and Numeric objects */
    /* ---------------------------------------------------------------------- */

    /* current front is empty */
    Work->fnpiv = 0 ;
    Work->fncols = 0 ;
    Work->fnrows = 0 ;
    Work->fnzeros = 0 ;

    Work->overlap = 0 ;
    Work->nz = nz ;
    Work->prior_element = EMPTY ;
    Work->ulen = 0 ;
    Work->llen = 0 ;
    Work->npiv = 0 ;
    Work->frontid = 0 ;

    Row_degree = Numeric->Rperm ;
    Col_degree = Numeric->Cperm ;
    /* Row_tuples = Numeric->Uip ; */
    Row_tlen   = Numeric->Uilen ;
    /* Col_tuples = Numeric->Lip ; */
    Col_tlen   = Numeric->Lilen ;

    Frpos = Work->Frpos ;
    Fcpos = Work->Fcpos ;
    Wp = Work->Wp ;

    /* D = Numeric->D ; */
    Upos = Numeric->Upos ;
    Lpos = Numeric->Lpos ;
    /* D [0..min(n_row,n_col)] = 0 ; no need to initialize this */

    for (row = 0 ; row <= n_row ; row++)
    {
	Lpos [row] = EMPTY ;
	/* Row_tuples [row] = 0 ; set in UMF_build_tuples */
	Row_degree [row] = 0 ;
	Row_tlen [row] = 0 ;
	/* Frpos [row] = EMPTY ;  do this later */
    }

    for (col = 0 ; col <= n_col ; col++)
    {
	Upos [col] = EMPTY ;
	/* Col_tuples [col] = 0 ; set in UMF_build_tuples */
	/* Col_degree [col] = 0 ; initialized below */
	Col_tlen [col] = 0 ;
	Fcpos [col] = EMPTY ;
    }

    for (i = 0 ; i <= nn ; i++)
    {
	Wp [i] = EMPTY ;
    }
    Work->Wpflag = -2 ;

    /* When cleared, Wp [0..nn] is < 0 and > wpflag. */
    /* In row search, Wp [col] is set to wpflag, which is negative. */
    /* In col search, Wp [row] is set to a position, which is >= 0. */
    /* Wp is cleared immediately after the row and col searches. */

    /* no need to initialize Wm, Wio, Woi, and Woo */

    /* clear the external degree counters */
    Work->cdeg0 = 1 ;
    Work->rdeg0 = 1 ;

    E = Work->E ;

    Numeric->n_row = n_row ;
    Numeric->n_col = n_col ;
    Numeric->npiv = 0 ;
    Numeric->nnzpiv = 0 ;
    Numeric->min_udiag = 0.0 ;
    Numeric->max_udiag = 0.0 ;
    Numeric->rcond = 0.0 ;
    Numeric->isize = 0 ;
    Numeric->nLentries = 0 ;
    Numeric->nUentries = 0 ;
    Numeric->lnz = 0 ;
    Numeric->unz = 0 ;
    Numeric->maxfrsize = 0 ;
    Numeric->flops = 0. ;

    /* ---------------------------------------------------------------------- */
    /* allocate the elements, copy the columns of A, and initialize degrees */
    /* ---------------------------------------------------------------------- */

    /* also apply the row and column pre-ordering.  */

    DEBUG3 (("LOAD_MATRIX:\n")) ;

    /* construct the inverse row permutation (use Frpos as temp) */
    for (newrow = 0 ; newrow < n_row ; newrow++)
    {
	oldrow = Rperm_init [newrow] ;
	ASSERT (oldrow >= 0 && oldrow < n_row) ;
	Frpos [oldrow] = newrow ;
    }

    e = 0 ;
    E [e] = 0 ;
    for (k = 0 ; k < n_col ; k++)
    {
	oldcol = Cperm_init [k] ;
	ASSERT (oldcol >= 0 && oldcol < n_col) ;
	p1 = Ap [oldcol] ;
	p2 = Ap [oldcol+1] ;
	deg = p2 - p1 ;
	Col_degree [k] = deg ;
	if (deg < 0)
	{
	    return (FALSE) ;		/* pattern changed */
	}
	else if (deg > 0)
	{
	    e++ ;
	    E [e] = UMF_mem_alloc_element (Numeric, deg, 1, &Rows, &Cols, &C,
	        &size, &ep) ;
	    if (E [e] <= 0)
	    {
	        /* We ran out of memory, which can only mean that */
	        /* the pattern (Ap and or Ai) has changed (gotten larger). */
	        DEBUG0 (("Pattern has gotten larger - alloc el. failed\n")) ;
	        DEBUG0 (("column k = "ID" size "ID"\n", k, size)) ;
	        return (FALSE) ;	/* pattern changed */
	    }
	    Cols [0] = k ;
	    ip = Rows ;
	    xp = C ;
	    ilast = -1 ;
	    for (p = p1 ; p < p2 ; p++)
	    {
		oldrow = Ai [p] ;
		if (oldrow <= ilast || oldrow >= n_row)
		{
		    return (FALSE) ;	/* pattern changed */
		}
		ilast = oldrow ;
		newrow = Frpos [oldrow] ;
		*ip++ = newrow ;
		ASSIGN (*xp, Ax [p], Az [p]) ;
		xp++ ;
		Row_degree [newrow]++ ;
	    }
	}
    }

    Work->nel = e ;
    Work->nelorig = e ;

    Col_degree [n_col] = 0 ;

    for (e++ ; e <= n_col + n_inner ; e++)
    {
	E [e] = 0 ;
    }

    /* Frpos no longer needed */
    for (row = 0 ; row <= n_row ; row++)
    {
	Frpos [row] = EMPTY ;
    }

    /* ---------------------------------------------------------------------- */
    /* build the tuple lists */
    /* ---------------------------------------------------------------------- */

    /* if the memory usage changes, then the pattern has changed */

    (void) UMF_tuple_lengths (Numeric, Work, &unused) ;
    if (!UMF_build_tuples (Numeric, Work))
    {
	/* We ran out of memory, which can only mean that */
	/* the pattern (Ap and or Ai) has changed (gotten larger). */
	DEBUG0 (("Pattern has gotten larger - build tuples failed\n")) ;
	return (FALSE) ;	/* pattern changed */
    }

    Numeric->init_usage = Numeric->max_usage ;
    if (Symbolic->num_mem_init_usage != Numeric->init_usage)
    {
	/* the pattern (Ap and or Ai) has changed */
	DEBUG0 (("Pattern has changed in size\n")) ;
	return (FALSE) ;	/* pattern changed */
    }

    /* NOTE:  this test does not detect all changes to Ap and/or Ai since the */
    /* last call to UMFPACK_*symbolic.  The pattern can change with the memory*/
    /* usage remaining the same. */

    /* ---------------------------------------------------------------------- */
    /* construct the row merge sets */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i <= Symbolic->nfr ; i++)
    {
	Work->Front_new1strow [i] = Symbolic->Front_1strow [i] ;
    }

#ifndef NDEBUG
    UMF_dump_rowmerge (Numeric, Symbolic, Work) ;
    DEBUG6 (("Column form of original matrix:\n")) ;
    UMF_dump_col_matrix (Ax,
#ifdef COMPLEX
    	Az,
#endif
	Ai, Ap, n_row, n_col, nz) ;
    UMF_gprob = gprob_save ;
    UMF_dump_memory (Numeric) ;
    UMF_dump_matrix (Numeric, Work, FALSE) ;
    DEBUG0 (("kernel init done...\n")) ;
#endif

    return (TRUE) ;

}

