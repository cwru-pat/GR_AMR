/* ========================================================================== */
/* === UMF_build_tuples ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Construct the tuple lists from a set of packed elements (no holes in
    elements, no internal or external fragmentation, and a packed (0..Work->nel)
    element name space).  Assume no tuple lists are currently allocated, but
    that the tuple lengths have been initialized by UMF_tuple_lengths.

    Returns TRUE if successful, FALSE if not enough memory.
*/

#include "umf_internal.h"
#include "umf_mem_alloc_tail_block.h"

GLOBAL Int UMF_build_tuples
(
    NumericType *Numeric,
    WorkType *Work
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int e, nrows, ncols, nel, *Rows, *Cols, row, col, n_row, n_col, *E,
	*Row_tuples, *Row_degree, *Row_tlen, tlen, tsize,
	*Col_tuples, *Col_degree, *Col_tlen ;
    Element *ep ;
    Unit *p ;
    Tuple tuple, *tp ;

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    E = Work->E ;
    Col_degree = Numeric->Cperm ;	/* for NON_PIVOTAL_COL macro */
    Row_degree = Numeric->Rperm ;	/* for NON_PIVOTAL_ROW macro */
    Row_tuples = Numeric->Uip ;
    Row_tlen   = Numeric->Uilen ;
    Col_tuples = Numeric->Lip ;
    Col_tlen   = Numeric->Lilen ;
    n_row = Work->n_row ;
    n_col = Work->n_col ;
    nel = Work->nel ;

    DEBUG3 (("BUILD_TUPLES: n_row "ID" n_col "ID" nel "ID"\n",
	n_row, n_col, nel)) ;

    /* ---------------------------------------------------------------------- */
    /* allocate space for the tuple lists */
    /* ---------------------------------------------------------------------- */

    /* Garbage collection and memory reallocation have already attempted to */
    /* ensure that there is enough memory for all the tuple lists.  If */
    /* memory allocation fails here, then there is nothing more to be done. */

    for (row = 0 ; row < n_row ; row++)
    {
	if (NON_PIVOTAL_ROW (row))
	{
	    tlen = MAX (4, Row_tlen [row] + 1) ;
	    tsize = UNITS (Tuple, tlen) ;
	    Row_tuples [row] = UMF_mem_alloc_tail_block (Numeric, tsize) ;
	    if (!Row_tuples [row])
	    {
		return (FALSE) ;	/* out of memory for row tuples */
	    }
	    Row_tlen [row] = 0 ;
	}
    }

    /* push on stack in reverse order, so column tuples are in the order */
    /* that they will be deleted. */
    for (col = n_col-1 ; col >= 0 ; col--)
    {
	if (NON_PIVOTAL_COL (col))
	{
	    tlen = MAX (4, Col_tlen [col] + 1) ;
	    tsize = UNITS (Tuple, tlen) ;
	    Col_tuples [col] = UMF_mem_alloc_tail_block (Numeric, tsize) ;
	    if (!Col_tuples [col])
	    {
		return (FALSE) ;	/* out of memory for col tuples */
	    }
	    Col_tlen [col] = 0 ;
	}
    }

#ifndef NDEBUG
    UMF_dump_memory (Numeric) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* create the tuple lists (exclude element 0) */
    /* ---------------------------------------------------------------------- */

    /* for all elements, in order of creation */
    for (e = 1 ; e <= nel ; e++)
    {
	DEBUG9 (("Adding tuples for element: "ID" at "ID"\n", e, E [e])) ;
	ASSERT (E [e]) ;	/* no external fragmentation */
	p = Numeric->Memory + E [e] ;
	GET_ELEMENT_PATTERN (ep, p, Cols, Rows, ncols) ;
	nrows = ep->nrows ;
	ASSERT (e != 0) ;
	ASSERT (e == 0 || nrows == ep->nrowsleft) ;
	ASSERT (e == 0 || ncols == ep->ncolsleft) ;
	tuple.e = e ;
	for (tuple.f = 0 ; tuple.f < ncols ; tuple.f++)
	{
	    col = Cols [tuple.f] ;
	    ASSERT (col >= 0 && col < n_col) ;
	    ASSERT (NON_PIVOTAL_COL (col)) ;
	    ASSERT (Col_tuples [col]) ;
	    tp = ((Tuple *) (Numeric->Memory + Col_tuples [col]))
		+ Col_tlen [col]++ ;
	    *tp = tuple ;
#ifndef NDEBUG
	    UMF_dump_rowcol (1, Numeric, Work, col, FALSE) ;
#endif
	}
	for (tuple.f = 0 ; tuple.f < nrows ; tuple.f++)
	{
	    row = Rows [tuple.f] ;
	    ASSERT (row >= 0 && row < n_row) ;
	    ASSERT (NON_PIVOTAL_COL (col)) ;
	    ASSERT (Row_tuples [row]) ;
	    tp = ((Tuple *) (Numeric->Memory + Row_tuples [row]))
		+ Row_tlen [row]++ ;
	    *tp = tuple ;
#ifndef NDEBUG
	    UMF_dump_rowcol (0, Numeric, Work, row, FALSE) ;
#endif
	}
    }

    /* ---------------------------------------------------------------------- */
    /* the tuple lists are now valid, and can be scanned */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    UMF_dump_memory (Numeric) ;
    UMF_dump_matrix (Numeric, Work, FALSE) ;
#endif
    DEBUG3 (("BUILD_TUPLES: done\n")) ;
    return (TRUE) ;
}

