/* ========================================================================== */
/* === UMF_garbage_collection =============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Compress the elements at the tail of Numeric->Memory, and delete the tuples.
    Elements are renumbered.  The new numbering space is compressed, and
    in the order of element creation (original elements of A first, followed
    by the new elements in the order that they were formed).
*/

#include "umf_internal.h"

GLOBAL void UMF_garbage_collection
(
    NumericType *Numeric,
    WorkType *Work
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int size, e, n_row, n_col, nrows, ncols, nrowsleft, ncolsleft, prevsize,
	csize, size2, i2, j2, i, j, cdeg, rdeg, *E, row, col, n_inner,
	*Rows, *Cols, *Rows2, *Cols2, nel, e2, *Row_tuples, *Col_tuples,
	*Row_degree, *Col_degree ;
    Entry *C, *C1, *C3, *C2 ;
    Unit *psrc, *pdest, *p, *pnext ;
    Element *epsrc, *epdest ;

#ifndef NDEBUG
    Int nmark ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    Col_degree = Numeric->Cperm ;	/* for NON_PIVOTAL_COL macro */
    Row_degree = Numeric->Rperm ;	/* for NON_PIVOTAL_ROW macro */
    Row_tuples = Numeric->Uip ;
    Col_tuples = Numeric->Lip ;
    E = Work->E ;
    n_row = Work->n_row ;
    n_col = Work->n_col ;
    n_inner = MIN (n_row, n_col) ;

    /* note that the tuple lengths (Col_tlen and Row_tlen) are updated, but */
    /* the tuple lists themselves are stale and are about to be destroyed */
    /* and recreated.  Do not attempt to scan them until they are recreated. */

#ifndef NDEBUG
    DEBUG0 (("::::GARBAGE COLLECTION::::\n")) ;
    UMF_dump_memory (Numeric) ;
#endif

    Numeric->ngarbage++ ;

    /* ---------------------------------------------------------------------- */
    /* delete the tuple lists by marking the blocks as free */
    /* ---------------------------------------------------------------------- */

    /* do not modify Row_tlen and Col_tlen */
    /* those are needed for UMF_build_tuples */

    for (row = 0 ; row < n_row ; row++)
    {
	if (NON_PIVOTAL_ROW (row) && Row_tuples [row])
	{
	    DEBUG2 (("row "ID" tuples "ID"\n", row, Row_tuples [row])) ;
	    p = Numeric->Memory + Row_tuples [row] - 1 ;
	    DEBUG2 (("Freeing tuple list row "ID", p-S "ID", size "ID"\n",
		row, p-Numeric->Memory, p->header.size)) ;
	    ASSERT (p->header.size > 0) ;
	    ASSERT (p >= Numeric->Memory + Numeric->itail) ;
	    ASSERT (p < Numeric->Memory + Numeric->size) ;
	    p->header.size = -p->header.size ;
	    Row_tuples [row] = 0 ;
	}
    }

    for (col = 0 ; col < n_col ; col++)
    {
	if (NON_PIVOTAL_COL (col) && Col_tuples [col])
	{
	    DEBUG2 (("col "ID" tuples "ID"\n", col, Col_tuples [col])) ;
	    p = Numeric->Memory + Col_tuples [col] - 1 ;
	    DEBUG2 (("Freeing tuple list col "ID", p-S "ID", size "ID"\n",
		col, p-Numeric->Memory, p->header.size)) ;
	    ASSERT (p->header.size > 0) ;
	    ASSERT (p >= Numeric->Memory + Numeric->itail) ;
	    ASSERT (p < Numeric->Memory + Numeric->size) ;
	    p->header.size = -p->header.size ;
	    Col_tuples [col] = 0 ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* mark the elements, and compress the name space */
    /* ---------------------------------------------------------------------- */

    nel = Work->nel ;

#ifndef NDEBUG
    nmark = 0 ;
#endif

    e2 = 0 ;
    ASSERT (!E [0]) ;	/* current version never uses element zero */

    for (e = 0 ; e <= nel ; e++) /* for all elements in order of creation */
    {
	if (E [e])
	{
	    psrc = Numeric->Memory + E [e] ;
	    psrc-- ;		/* get the header of this block */
	    if (e > 0)
	    {
		e2++ ;	/* do not renumber element zero */
	    }
	    ASSERT (psrc->header.size > 0) ;
	    psrc->header.size = e2  ;	/* store the new name in the header */
#ifndef NDEBUG
	    nmark++ ;
#endif
	    DEBUG7 ((ID":: Mark e "ID" at psrc-S "ID", new e "ID"\n",
		nmark, e, psrc-Numeric->Memory, e2)) ;
	    E [e] = 0 ;

	    if (e == Work->nelorig)
	    {
		/* this is the highest numbered original element */
		Work->nelorig = e2 ;
	    }
	    if (e == Work->prior_element)
	    {
		Work->prior_element = e2 ;
	    }
	}
    }

    /* all 1..e2 are now in use (element zero may or may not be in use) */
    Work->nel = e2 ;
    nel = Work->nel ;

#ifndef NDEBUG
    for (e = 0 ; e <= n_col + n_inner ; e++)
    {
	ASSERT (!E [e]) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* compress the elements */
    /* ---------------------------------------------------------------------- */

    /* point to tail marker block of size 1 + header */
    psrc = Numeric->Memory + Numeric->size - 2 ;
    pdest = psrc ;
    prevsize = psrc->header.prevsize ;
    DEBUG7 (("Starting the compression:\n")) ;

    while (prevsize > 0)
    {

	/* ------------------------------------------------------------------ */
	/* move up to the next element above the current header, and */
	/* get the element name and size */
	/* (if it is an element, the name will be positive) */
	/* ------------------------------------------------------------------ */

	size = prevsize ;
	psrc -= (size + 1) ;
	e = psrc->header.size ;
	prevsize = psrc->header.prevsize ;
	/* top block at tail has prevsize of 0 */

	/* a free block will have a negative size, so skip it */
	/* otherwise, if size > 0, it holds the element name, not the size */

	DEBUG8 (("psrc-S: "ID" prevsize: "ID" size: "ID, psrc-Numeric->Memory,
	prevsize, size)) ;

	if (e >= 0)
	{

	    /* -------------------------------------------------------------- */
	    /* this is an element, compress and move from psrc down to pdest */
	    /* -------------------------------------------------------------- */

#ifndef NDEBUG
	    DEBUG7 (("\n")) ;
	    DEBUG7 ((ID":: Move element "ID": from: "ID" \n",
		nmark, e, psrc-Numeric->Memory)) ;
	    nmark-- ;
	    ASSERT (e <= nel) ;
	    ASSERT (E [e] == 0) ;
#endif

	    /* -------------------------------------------------------------- */
	    /* get the element scalars, and pointers to C, Rows, and Cols: */
	    /* -------------------------------------------------------------- */

	    p = psrc + 1 ;
	    GET_ELEMENT (epsrc, p, Cols, Rows, ncols, nrows, C) ;
	    nrowsleft = epsrc->nrowsleft ;
	    ncolsleft = epsrc->ncolsleft ;
	    cdeg = epsrc->cdeg ;
	    rdeg = epsrc->rdeg ;

#ifndef NDEBUG
	    DEBUG7 ((" nrows "ID" nrowsleft "ID"\n", nrows, nrowsleft)) ;
	    DEBUG7 ((" ncols "ID" ncolsleft "ID"\n", ncols, ncolsleft)) ;
	    DEBUG8 ((" Rows:")) ;
	    for (i = 0 ; i < nrows ; i++) DEBUG8 ((" "ID, Rows [i])) ;
	    DEBUG8 (("\n Cols:")) ;
	    for (j = 0 ; j < ncols ; j++) DEBUG8 ((" "ID, Cols [j])) ;
	    DEBUG8 (("\n")) ;
#endif

	    /* -------------------------------------------------------------- */
	    /* determine the layout of the new element */
	    /* -------------------------------------------------------------- */

	    csize = nrowsleft * ncolsleft ;
	    size2 = UNITS (Element, 1)
		  + UNITS (Int, nrowsleft + ncolsleft)
		  + UNITS (Entry, csize) ;

	    DEBUG7 (("Old size "ID" New size "ID"\n", size, size2)) ;

	    pnext = pdest ;
	    pnext->header.prevsize = size2 ;
	    pdest -= (size2 + 1) ;
	    DEBUG7 (("Move element from psrc-S "ID" to pdest-S "ID"\n",
		psrc-Numeric->Memory, pdest-Numeric->Memory));

	    ASSERT (e > 0 || size == size2) ;
	    ASSERT (size2 <= size) ;
	    ASSERT (psrc + 1 + size <= pnext) ;
	    ASSERT (psrc <= pdest) ;

	    p = pdest + 1 ;
	    epdest = (Element *) p ;
	    p += UNITS (Element, 1) ;
	    Cols2 = (Int *) p ;
	    Rows2 = Cols2 + ncolsleft ;
	    p += UNITS (Int, nrowsleft + ncolsleft) ;
	    C2 = (Entry *) p ;

	    ASSERT (epdest >= epsrc) ;
	    ASSERT (Rows2 >= Rows) ;
	    ASSERT (Cols2 >= Cols) ;
	    ASSERT (C2 >= C) ;
	    ASSERT (p + UNITS (Entry, csize) == pnext) ;

	    /* -------------------------------------------------------------- */
	    /* move the contribution block */
	    /* -------------------------------------------------------------- */

	    /* overlap = psrc + size + 1 > pdest ; */

	    if (nrowsleft < nrows || ncolsleft < ncols)
	    {

		/* ---------------------------------------------------------- */
		/* compress contribution block in place prior to moving it */
		/* ---------------------------------------------------------- */

		DEBUG7 (("Compress C in place prior to move (incl unused):\n"));
#ifndef NDEBUG
		UMF_dump_dense (C, nrows, nrows, ncols) ;
#endif
		C1 = C ;
		C3 = C ;
		for (j = 0 ; j < ncols ; j++)
		{
		    if (Cols [j] >= 0)
		    {
			for (i = 0 ; i < nrows ; i++)
			{
			    if (Rows [i] >= 0)
			    {
				*C3++ = C1 [i] ;
			    }
			}
		    }
		    C1 += nrows ;
		}
		ASSERT (C3-C == csize) ;
		DEBUG8 (("Newly compressed contrib. block (all in use):\n")) ;
#ifndef NDEBUG
		UMF_dump_dense (C, nrowsleft, nrowsleft, ncolsleft) ;
#endif
	    }

	    /* shift the contribution block down */
	    C += csize ;
	    C2 += csize ;
	    for (i = 0 ; i < csize ; i++)
	    {
		*--C2 = *--C ;
	    }

	    /* -------------------------------------------------------------- */
	    /* move the row indices */
	    /* -------------------------------------------------------------- */

	    i2 = nrowsleft ;
	    for (i = nrows - 1 ; i >= 0 ; i--)
	    {
		ASSERT (Rows2+i2 >= Rows+i) ;
		if (Rows [i] >= 0)
		{
		    Rows2 [--i2] = Rows [i] ;
		}
	    }
	    ASSERT (i2 == 0) ;

	    j2 = ncolsleft ;
	    for (j = ncols - 1 ; j >= 0 ; j--)
	    {
		ASSERT (Cols2+j2 >= Cols+j) ;
		if (Cols [j] >= 0)
		{
		    Cols2 [--j2] = Cols [j] ;
		}
	    }
	    ASSERT (j2 == 0) ;

	    /* -------------------------------------------------------------- */
	    /* construct the new header */
	    /* -------------------------------------------------------------- */

	    /* E [0...e] is now valid */
	    E [e] = (pdest + 1) - Numeric->Memory ;
	    epdest = (Element *) (pdest + 1) ;

	    epdest->next = EMPTY ;	/* destroys the son list */
	    epdest->ncols = ncolsleft ;
	    epdest->nrows = nrowsleft ;
	    epdest->ncolsleft = ncolsleft ;
	    epdest->nrowsleft = nrowsleft ;
	    epdest->rdeg = rdeg ;
	    epdest->cdeg = cdeg ;

	    ASSERT (size2 <= size) ;
	    pdest->header.prevsize = 0 ;
	    pdest->header.size = size2 ;

	    DEBUG7 (("After moving it:\n")) ;
#ifndef NDEBUG
	    UMF_dump_element (Numeric, Work, e, FALSE) ;
#endif
	}

#ifndef NDEBUG
	else
	{
	    DEBUG8 ((" free\n")) ;
	}
#endif
	DEBUG7 (("psrc "ID"  tail "ID"\n",
	psrc-Numeric->Memory, Numeric->itail)) ;
    }

    ASSERT (psrc == Numeric->Memory + Numeric->itail) ;
    ASSERT (nmark == 0) ;

    /* ---------------------------------------------------------------------- */
    /* final tail pointer */
    /* ---------------------------------------------------------------------- */

    ASSERT (pdest >= Numeric->Memory + Numeric->itail) ;
    Numeric->itail = pdest - Numeric->Memory ;
    pdest->header.prevsize = 0 ;
    Numeric->ibig = EMPTY ;
    Numeric->tail_usage = Numeric->size - Numeric->itail ;

    /* ---------------------------------------------------------------------- */
    /* clear the unused E [nel+1 .. n_col + n_inner] */
    /* ---------------------------------------------------------------------- */

    for (e = nel+1 ; e <= n_col + n_inner ; e++)
    {
	E [e] = 0 ;
    }

#ifndef NDEBUG
    UMF_dump_packed_memory (Numeric, Work) ;
#endif

    DEBUG8 (("::::GARBAGE COLLECTION DONE::::\n")) ;
}

