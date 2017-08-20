/* ========================================================================== */
/* === UMF_create_element =================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Factorization of a frontal matrix is complete.  Create a new element for
    later assembly into a subsequent frontal matrix.  Returns TRUE if
    successful, FALSE if out of memory.
*/

#include "umf_internal.h"
#include "umf_mem_alloc_element.h"
#include "umf_mem_alloc_tail_block.h"
#include "umf_mem_free_tail_block.h"
#include "umf_get_memory.h"
#include "umf_blas3_update.h"

GLOBAL Int UMF_create_element
(
    NumericType *Numeric,
    WorkType *Work
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int j, col, row, *Fcols, *Frows, fnrows, fncols, *Cols, len, needunits, t1,
	t2, size, e, i, *E, *Fcpos, *Frpos, *Rows, eloc, fnrows_max, f,
	got_memory, *Row_tuples, *Row_degree, *Row_tlen, *Col_tuples, max_mark,
	*Col_degree, *Col_tlen, nn, n_row, n_col ;
    Entry *C, *Fcol ;
    Element *ep ;
    Unit *p, *Memory ;
    Tuple *tp, *tp1, *tp2, tuple, *tpend ;
#ifndef NDEBUG
    double gprob_save ;
#endif

    /* ---------------------------------------------------------------------- */
    /* apply pending updates, if any */
    /* ---------------------------------------------------------------------- */

    if (Work->fnpiv > 0)
    {
	UMF_blas3_update (Work) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    DEBUG2 (("FRONTAL WRAPUP\n")) ;
    Col_degree = Numeric->Cperm ;
    Row_degree = Numeric->Rperm ;
    Row_tuples = Numeric->Uip ;
    Row_tlen   = Numeric->Uilen ;
    Col_tuples = Numeric->Lip ;
    Col_tlen   = Numeric->Lilen ;
    n_row = Work->n_row ;
    n_col = Work->n_col ;
    nn = MAX (n_row, n_col) ;
    Fcols = Work->Fcols ;
    Frows = Work->Frows ;
    Fcpos = Work->Fcpos ;
    Frpos = Work->Frpos ;
    Memory = Numeric->Memory ;
    fncols = Work->fncols ;
    fnrows = Work->fnrows ;
    fnrows_max = Work->fnrows_max ;

    tp = (Tuple *) NULL ;
    tp1 = (Tuple *) NULL ;
    tp2 = (Tuple *) NULL ;

    /* ---------------------------------------------------------------------- */
    /* add the current frontal matrix to the degrees of each column, */
    /* and place the modified columns back in the degree lists */
    /* ---------------------------------------------------------------------- */

    for (j = 0 ; j < fncols ; j++)
    {
	/* add the current frontal matrix to the degree */
	ASSERT (Fcols [j] >= 0 && Fcols [j] < n_col) ;
	Col_degree [Fcols [j]] += fnrows ;
    }

    /* ---------------------------------------------------------------------- */
    /* add the current frontal matrix to the degrees of each row */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < fnrows ; i++)
    {
	/* add the current frontal matrix to the degree */
	ASSERT (Frows [i] >= 0 && Frows [i] < n_row) ;
	Row_degree [Frows [i]] += fncols ;
    }

    /* ---------------------------------------------------------------------- */
    /* Reset the external degree counters */
    /* ---------------------------------------------------------------------- */

    E = Work->E ;
    max_mark = MAX_MARK (nn) ;

    if (!Work->pivcol_in_front)
    {
	/* clear the external column degrees. no more Usons of current front */
	Work->cdeg0 += (nn + 1) ;
	if (Work->cdeg0 >= max_mark)
	{
	    /* guard against integer overflow.  This is very rare */
	    Work->cdeg0 = 1 ;
	    for (e = 0 ; e <= Work->nel ; e++)
	    {
		if (E [e])
		{
		    ep = (Element *) (Memory + E [e]) ;
		    ep->cdeg = 0 ;
		}
	    }
	}
    }

    if (!Work->pivrow_in_front)
    {
	/* clear the external row degrees.  no more Lsons of current front */
	Work->rdeg0 += (nn + 1) ;
	if (Work->rdeg0 >= max_mark)
	{
	    /* guard against integer overflow.  This is very rare */
	    Work->rdeg0 = 1 ;
	    for (e = 0 ; e <= Work->nel ; e++)
	    {
		if (E [e])
		{
		    ep = (Element *) (Memory + E [e]) ;
		    ep->rdeg = 0 ;
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* clear row/col offsets */
    /* ---------------------------------------------------------------------- */

    if (!Work->pivrow_in_front)
    {
	for (j = 0 ; j < fncols ; j++)
	{
	    Fcpos [Fcols [j]] = EMPTY ;
	}
    }

    if (!Work->pivcol_in_front)
    {
	for (i = 0 ; i < fnrows ; i++)
	{
	    Frpos [Frows [i]] = EMPTY ;
	}
    }

    if (fncols <= 0 || fnrows <= 0)
    {
	/* no element to create */
	DEBUG2 (("Element evaporation\n")) ;
	Work->prior_element = EMPTY ;
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* create element for later assembly */
    /* ---------------------------------------------------------------------- */

    needunits = 0 ;
    got_memory = FALSE ;
    eloc = UMF_mem_alloc_element (Numeric, fnrows, fncols, &Rows, &Cols, &C,
	&needunits, &ep) ;

    if (!eloc)
    {
	/* do garbage collection, realloc, and try again */
	if (!UMF_get_memory (Numeric, Work, needunits))
	{
	    return (FALSE) ;	/* out of memory */
	}
	got_memory = TRUE ;
	Memory = Numeric->Memory ;
#ifndef NDEBUG
	gprob_save = UMF_gprob ;
	UMF_gprob = 0 ;
#endif
	eloc = UMF_mem_alloc_element (Numeric, fnrows, fncols, &Rows, &Cols, &C,
	    &needunits, &ep) ;
	ASSERT (eloc >= 0) ;
	if (!eloc)
	{
	    return (FALSE) ;	/* out of memory */
	}
#ifndef NDEBUG
	UMF_gprob = gprob_save ;
#endif
    }

    e = ++(Work->nel) ;	/* get the name of this new frontal matrix */
    Work->prior_element = e ;
    DEBUG8 (("wrapup e "ID" nel "ID"\n", e, Work->nel)) ;

    ASSERT (e > 0 && e <= n_col + MIN (n_row,n_col)) ;
    ASSERT (E [e] == 0) ;
    E [e] = eloc ;

    if (Work->pivcol_in_front)
    {
	/* the new element is a Uson of the next frontal matrix */
	ep->cdeg = Work->cdeg0 ;
    }

    if (Work->pivrow_in_front)
    {
	/* the new element is an Lson of the next frontal matrix */
	ep->rdeg = Work->rdeg0 ;
    }

    /* ---------------------------------------------------------------------- */
    /* copy frontal matrix into the new element */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < fnrows ; i++)
    {
	Rows [i] = Frows [i] ;
    }
    for (i = 0 ; i < fncols ; i++)
    {
	Cols [i] = Fcols [i] ;
    }
    Fcol = Work->Fx ;
    for (j = 0 ; j < fncols ; j++)
    {
	for (i = 0 ; i < fnrows ; i++)
	{
	    C [i] = Fcol [i] ;
	}
	Fcol += fnrows_max ;
	C += fnrows ;
    }

    /* ---------------------------------------------------------------------- */
    /* add tuples for the new element */
    /* ---------------------------------------------------------------------- */

    tuple.e = e ;

    if (got_memory)
    {

	/* ------------------------------------------------------------------ */
	/* UMF_get_memory ensures enough space exists for each new tuple */
	/* ------------------------------------------------------------------ */

	/* place (e,f) in the element list of each column */
	for (tuple.f = 0 ; tuple.f < fncols ; tuple.f++)
	{
	    col = Fcols [tuple.f] ;
	    ASSERT (col >= 0 && col < n_col) ;
	    ASSERT (NON_PIVOTAL_COL (col)) ;
	    ASSERT (Col_tuples [col]) ;
	    tp = ((Tuple *) (Memory + Col_tuples [col])) + Col_tlen [col]++ ;
	    *tp = tuple ;
	}

	/* ------------------------------------------------------------------ */

	/* place (e,f) in the element list of each row */
	for (tuple.f = 0 ; tuple.f < fnrows ; tuple.f++)
	{
	    row = Frows [tuple.f] ;
	    ASSERT (row >= 0 && row < n_row) ;
	    ASSERT (NON_PIVOTAL_ROW (row)) ;
	    ASSERT (Row_tuples [row]) ;
	    tp = ((Tuple *) (Memory + Row_tuples [row])) + Row_tlen [row]++ ;
	    *tp = tuple ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* might not have enough space for each tuple */
	/* ------------------------------------------------------------------ */

	/* ------------------------------------------------------------------ */
	/* place (e,f) in the element list of each column */
	/* ------------------------------------------------------------------ */

	for (tuple.f = 0 ; tuple.f < fncols ; tuple.f++)
	{
	    col = Fcols [tuple.f] ;
	    ASSERT (col >= 0 && col < n_col) ;
	    ASSERT (NON_PIVOTAL_COL (col)) ;
	    t1 = Col_tuples [col] ;
	    DEBUG1 (("Placing on col:"ID" , tuples at "ID"\n",
		col, Col_tuples [col])) ;

	    size = 0 ;
	    len = 0 ;

	    if (t1)
	    {
		p = Memory + t1 ;
		tp = (Tuple *) p ;
		size = GET_BLOCK_SIZE (p) ;
		len = Col_tlen [col] ;
		tp2 = tp + len ;
	    }

	    needunits = UNITS (Tuple, len + 1) ;
	    DEBUG1 (("len: "ID" size: "ID" needunits: "ID"\n",
		len, size, needunits));

	    if (needunits > size && t1)
	    {
		/* prune the tuples */
		tp1 = tp ;
		tp2 = tp ;
		tpend = tp + len ;
		for ( ; tp < tpend ; tp++)
		{
		    e = tp->e ;
		    ASSERT (e > 0 && e <= Work->nel) ;
		    if (!E [e]) continue ;   /* element already deallocated */
		    f = tp->f ;
		    p = Memory + E [e] ;
		    ep = (Element *) p ;
		    p += UNITS (Element, 1) ;
		    Cols = (Int *) p ;
		    ;
		    if (Cols [f] == EMPTY) continue ;	/* already assembled */
		    ASSERT (col == Cols [f]) ;
		    *tp2++ = *tp ;	/* leave the tuple in the list */
		}
		len = tp2 - tp1 ;
		Col_tlen [col] = len ;
		needunits = UNITS (Tuple, len + 1) ;
	    }

	    if (needunits > size)
	    {
		/* no room exists - reallocate elsewhere */
		DEBUG1 (("REALLOCATE Col: "ID", size "ID" to "ID"\n",
		    col, size, 2*needunits)) ;
#ifndef NDEBUG
		UMF_allocfail = FALSE ;
		if (UMF_gprob > 0)  /* a double relop, but ignore NaN case */
		{
		    double rrr = ((double) (rand ( ))) /
		    	(((double) RAND_MAX) + 1) ;
		    DEBUG1 (("Check random %e %e\n", rrr, UMF_gprob)) ;
		    UMF_allocfail = rrr < UMF_gprob ;
		    if (UMF_allocfail) DEBUG1 (("Random gar. (col tuple)\n")) ;
		}
#endif
		needunits = MIN (2*needunits, (Int) UNITS (Tuple, nn)) ;
		t2 = UMF_mem_alloc_tail_block (Numeric, needunits) ;
		if (!t2)
		{
		    /* get memory, reconstruct all tuple lists, and return */
		    return (UMF_get_memory (Numeric, Work, 0)) ;
		}
		Col_tuples [col] = t2 ;
		tp2 = (Tuple *) (Memory + t2) ;
		if (t1)
		{
		    for (i = 0 ; i < len ; i++)
		    {
			*tp2++ = *tp1++ ;
		    }
		    UMF_mem_free_tail_block (Numeric, t1) ;
		}
	    }

	    /* place the new (e,f) tuple in the element list of the column */
	    Col_tlen [col]++ ;
	    *tp2 = tuple ;
	}

	/* ------------------------------------------------------------------ */
	/* place (e,f) in the element list of each row */
	/* ------------------------------------------------------------------ */

	for (tuple.f = 0 ; tuple.f < fnrows ; tuple.f++)
	{
	    row = Frows [tuple.f] ;
	    ASSERT (row >= 0 && row < n_row) ;
	    ASSERT (NON_PIVOTAL_ROW (row)) ;
	    t1 = Row_tuples [row] ;
	    DEBUG1 (("Placing on row:"ID" , tuples at "ID"\n",
		row, Row_tuples [row])) ;

	    size = 0 ;
	    len = 0 ;
	    if (t1)
	    {
		p = Memory + t1 ;
		tp = (Tuple *) p ;
		size = GET_BLOCK_SIZE (p) ;
		len = Row_tlen [row] ;
		tp2 = tp + len ;
	    }

	    needunits = UNITS (Tuple, len + 1) ;
	    DEBUG1 (("len: "ID" size: "ID" needunits: "ID"\n",
		len, size, needunits)) ;

	    if (needunits > size && t1)
	    {
		/* prune the tuples */
		tp1 = tp ;
		tp2 = tp ;
		tpend = tp + len ;
		for ( ; tp < tpend ; tp++)
		{
		    e = tp->e ;
		    ASSERT (e > 0 && e <= Work->nel) ;
		    if (!E [e])
		    {
			continue ;	/* element already deallocated */
		    }
		    f = tp->f ;
		    p = Memory + E [e] ;
		    ep = (Element *) p ;
		    p += UNITS (Element, 1) ;
		    Cols = (Int *) p ;
		    Rows = Cols + (ep->ncols) ;
		    if (Rows [f] == EMPTY) continue ;	/* already assembled */
		    ASSERT (row == Rows [f]) ;
		    *tp2++ = *tp ;	/* leave the tuple in the list */
		}
		len = tp2 - tp1 ;
		Row_tlen [row] = len ;
		needunits = UNITS (Tuple, len + 1) ;
	    }

	    if (needunits > size)
	    {
		/* no room exists - reallocate elsewhere */
		DEBUG1 (("REALLOCATE Row: "ID", size "ID" to "ID"\n",
		    row, size, 2*needunits)) ;
#ifndef NDEBUG
		UMF_allocfail = FALSE ;
		if (UMF_gprob > 0)  /* a double relop, but ignore NaN case */
		{
		    double rrr = ((double) (rand ( ))) /
		    	(((double) RAND_MAX) + 1) ;
		    DEBUG1 (("Check random %e %e\n", rrr, UMF_gprob)) ;
		    UMF_allocfail = rrr < UMF_gprob ;
		    if (UMF_allocfail) DEBUG1 (("Random gar. (row tuple)\n")) ;
		}
#endif
		needunits = MIN (2*needunits, (Int) UNITS (Tuple, nn)) ;
		t2 = UMF_mem_alloc_tail_block (Numeric, needunits) ;
		if (!t2)
		{
		    /* get memory, reconstruct all tuple lists, and return */
		    return (UMF_get_memory (Numeric, Work, 0)) ;
		}
		Row_tuples [row] = t2 ;
		tp2 = (Tuple *) (Memory + t2) ;
		if (t1)
		{
		    for (i = 0 ; i < len ; i++)
		    {
			*tp2++ = *tp1++ ;
		    }
		    UMF_mem_free_tail_block (Numeric, t1) ;
		}
	    }

	    /* place the new (e,f) tuple in the element list of the row */
	    Row_tlen [row]++ ;
	    *tp2 = tuple ;
	}

    }

    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    DEBUG1 (("Done extending\nFINAL: element row pattern: len="ID"\n", fncols));
    for (j = 0 ; j < fncols ; j++) DEBUG1 ((""ID"\n", Fcols [j])) ;
    DEBUG1 (("FINAL: element col pattern:  len="ID"\n", fnrows)) ;
    for (j = 0 ; j < fnrows ; j++) DEBUG1 ((""ID"\n", Frows [j])) ;
    for (j = 0 ; j < fncols ; j++)
    {
	col = Fcols [j] ;
	ASSERT (col >= 0 && col < n_col) ;
	UMF_dump_rowcol (1, Numeric, Work, col, TRUE) ;
    }
    for (j = 0 ; j < fnrows ; j++)
    {
	row = Frows [j] ;
	ASSERT (row >= 0 && row < n_row) ;
	UMF_dump_rowcol (0, Numeric, Work, row, TRUE) ;
    }

    if (n_row < 1000 && n_col < 1000)
    {
	UMF_dump_memory (Numeric) ;
    }

    DEBUG1 (("New element, after filling with stuff: "ID"\n", e)) ;
    UMF_dump_element (Numeric, Work, e, TRUE) ;

    if (nn < 1000)
    {
	DEBUG4 (("Matrix dump, after New element: "ID"\n", e)) ;
	UMF_dump_matrix (Numeric, Work, TRUE) ;
    }

    DEBUG3 (("FRONTAL WRAPUP DONE\n")) ;
#endif

    return (TRUE) ;
}

