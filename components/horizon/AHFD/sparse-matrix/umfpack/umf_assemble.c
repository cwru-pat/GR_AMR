/* ========================================================================== */
/* === UMF_assemble ========================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*  Degree update and numerical assembly */

#include "umf_internal.h"
#include "umf_mem_free_tail_block.h"

GLOBAL void UMF_assemble
(
    NumericType *Numeric,
    WorkType *Work
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int e, i, row, col, i2, nrows, ncols, f, tpi, extcdeg, extrdeg, rdeg0,
	cdeg0, son_list, next, scan_pivrow, scan_pivcol, nrows_to_assemble,
	ncols_to_assemble, ngetrows, j, j2,
	nrowsleft,	/* number of rows remaining in C */
	ncolsleft,	/* number of columns remaining in C  */
	prior_Lson, prior_Uson, *E, *Cols, *Rows, *Wm, *Woo,
	*Row_tuples, *Row_degree, *Row_tlen,
	*Col_tuples, *Col_degree, *Col_tlen ;
    Unit *Memory, *p ;
    Element *ep ;
    Tuple *tp, *tp1, *tp2, *tpend ;

    Entry
	*C,		/* a pointer into the contribution block */
	*Fx,		/* Fx [0..fnrows_max-1, 0..fncols_max-1] working array*/
	*Fcol,		/* a column of F */
	*Frow ;		/* a row of F */

    Int
	*Frows,		/* Frows [0.. ]: row indices of F */
	*Fcols,		/* Fcols [0.. ]: column indices of F */
	*Frpos,
	*Fcpos,
	fnrows,		/* number of rows in contribution block in F */
	fncols ;	/* number of columns in contribution block in F */

#ifndef NDEBUG
    Int j3, n_row, n_col ;
    n_row = Work->n_row ;
    n_col = Work->n_col ;
    DEBUG3 (("::Assemble SCANS 1-4\n")) ;
    UMF_dump_current_front (Numeric, Work, TRUE) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    fncols = Work->fncols ;
    fnrows = Work->fnrows ;
    Fcols = Work->Fcols ;
    Frows = Work->Frows ;
    Fcpos = Work->Fcpos ;
    Frpos = Work->Frpos ;
    Fx = Work->Fx ;
    Col_degree = Numeric->Cperm ;
    Row_degree = Numeric->Rperm ;
    Row_tuples = Numeric->Uip ;
    Row_tlen   = Numeric->Uilen ;
    Col_tuples = Numeric->Lip ;
    Col_tlen   = Numeric->Lilen ;
    E = Work->E ;
    Memory = Numeric->Memory ;
    Wm = Work->Wm ;
    Woo = Work->Woo ;		/* must be of size at least fnrows_max */
    rdeg0 = Work->rdeg0 ;
    cdeg0 = Work->cdeg0 ;

#ifndef NDEBUG
    DEBUG6 (("============================================\n")) ;
    DEBUG6 (("Degree update, assembly.\n")) ;
    DEBUG6 (("pivot row pattern:  fncols="ID"\n", fncols)) ;
    for (j3 = 0 ; j3 < fncols ; j3++) DEBUG6 ((ID, Fcols [j3])) ;
    DEBUG6 (("\npivot col pattern:  fnrows="ID"\n", fnrows)) ;
    for (j3 = 0 ; j3 < fnrows ; j3++) DEBUG6 ((ID, Frows [j3])) ;
    DEBUG6 (("\n")) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* determine the largest actual frontal matrix size */
    /* ---------------------------------------------------------------------- */

    Numeric->maxfrsize = MAX (Numeric->maxfrsize,
	(fnrows + Work->fnpiv) * (fncols + Work->fnpiv)) ;

    /* ---------------------------------------------------------------------- */
    /* assemble from prior elements into the current frontal matrix */
    /* ---------------------------------------------------------------------- */

    DEBUG2 (("New assemble start [\n")) ;

    /* Currently no rows or columns are marked.  No elements are scanned, */
    /* that is, (ep->next == EMPTY) is true for all elements */

    son_list = 0 ;	/* start creating son_list [ */

    /* ---------------------------------------------------------------------- */
    /* determine if most recent element is Lson or Uson of current front */
    /* ---------------------------------------------------------------------- */

    if (!Work->do_extend)
    {
	prior_Uson = ( Work->pivcol_in_front && !Work->pivrow_in_front) ;
	prior_Lson = (!Work->pivcol_in_front &&  Work->pivrow_in_front) ;
	if (prior_Uson || prior_Lson)
	{
	    e = Work->prior_element ;
	    if (e != EMPTY)
	    {
		ASSERT (E [e]) ;
		p = Memory + E [e] ;
		ep = (Element *) p ;
		ep->next = son_list ;
		son_list = e ;
#ifndef NDEBUG
		DEBUG2 (("e "ID" is Prior son "ID" "ID"\n",
		    e, prior_Uson, prior_Lson)) ;
		UMF_dump_element (Numeric, Work, e, FALSE) ;
#endif
		ASSERT (E [e]) ;
	    }
	}
    }
    Work->prior_element = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* SCAN1-row:  scan the element lists of each new row in the pivot col */
    /* and compute the external column degree for each frontal */
    /* ---------------------------------------------------------------------- */

    scan_pivrow = Work->scan_pivrow ;

    for (i2 = Work->fscan_row ; i2 < fnrows || scan_pivrow ; )
    {
	/* Get a row */
	if (scan_pivrow)
	{
	    /* pivot row is new to this front.  Scan it */
	    row = Work->pivrow ;
	    scan_pivrow = FALSE ;
	}
	else
	{
	    row = Frows [i2++] ;
	}

	DEBUG6 (("SCAN1-row: "ID"\n", row)) ;
#ifndef NDEBUG
	UMF_dump_rowcol (0, Numeric, Work, row, FALSE) ;
#endif

	ASSERT (NON_PIVOTAL_ROW (row)) ;
	tpi = Row_tuples [row] ;
	if (!tpi) continue ;
	tp = (Tuple *) (Memory + tpi) ;
	tp1 = tp ;
	tp2 = tp ;
	tpend = tp + Row_tlen [row] ;
	for ( ; tp < tpend ; tp++)
	{
	    e = tp->e ;
	    ASSERT (e > 0 && e <= Work->nel) ;
	    if (!E [e]) continue ;	/* element already deallocated */
	    f = tp->f ;
	    p = Memory + E [e] ;
	    ep = (Element *) p ;
	    p += UNITS (Element, 1) ;
	    Rows = ((Int *) p) + ep->ncols ;
	    if (Rows [f] == EMPTY) continue ;	/* row already assembled */
	    ASSERT (row == Rows [f]) ;

	    if (ep->cdeg < cdeg0)
	    {
		/* first time seen in scan1-row */
		ep->cdeg = ep->nrowsleft + cdeg0 ;
		DEBUG6 (("e "ID" First seen: cdeg: "ID" ", e, ep->cdeg-cdeg0)) ;
		ASSERT (ep->ncolsleft > 0 && ep->nrowsleft > 0) ;
	    }

	    ep->cdeg-- ;	/* decrement external column degree */
	    DEBUG6 (("e "ID" New ext col deg: "ID"\n", e, ep->cdeg - cdeg0)) ;

	    /* this element is not yet in the new son list */
	    if (ep->cdeg == cdeg0 && ep->next == EMPTY)
	    {
		/* A new LUson or Uson has been found */
		ep->next = son_list ;
		son_list = e ;
	    }

	    ASSERT (ep->cdeg >= cdeg0) ;
	    *tp2++ = *tp ;	/* leave the tuple in the list */
	}
	Row_tlen [row] = tp2 - tp1 ;
    }

    /* ---------------------------------------------------------------------- */
    /* SCAN1-col:  scan the element lists of each new col in the pivot row */
    /*	 and compute the external row degree for each frontal */
    /* ---------------------------------------------------------------------- */

    scan_pivcol = Work->scan_pivcol ;

    for (j2 = Work->fscan_col ; j2 < fncols || scan_pivcol ; )
    {
	/* Get a column */
	if (scan_pivcol)
	{
	    /* pivot col is new to this front.  Scan it */
	    col = Work->pivcol ;
	    scan_pivcol = FALSE ;
	}
	else
	{
	    col = Fcols [j2++] ;
	}
	ASSERT (col >= 0 && col < n_col) ;

	DEBUG6 (("SCAN 1-col: "ID"\n", col)) ;
#ifndef NDEBUG
	UMF_dump_rowcol (1, Numeric, Work, col, FALSE) ;
#endif

	ASSERT (NON_PIVOTAL_COL (col)) ;
	tpi = Col_tuples [col] ;
	if (!tpi) continue ;
	tp = (Tuple *) (Memory + tpi) ;
	tp1 = tp ;
	tp2 = tp ;
	tpend = tp + Col_tlen [col] ;
	for ( ; tp < tpend ; tp++)
	{
	    e = tp->e ;
	    ASSERT (e > 0 && e <= Work->nel) ;
	    if (!E [e]) continue ;	/* element already deallocated */
	    f = tp->f ;
	    p = Memory + E [e] ;
	    ep = (Element *) p ;
	    p += UNITS (Element, 1) ;
	    Cols = (Int *) p ;
	    if (Cols [f] == EMPTY) continue ;	/* column already assembled */
	    ASSERT (col == Cols [f]) ;

	    if (ep->rdeg < rdeg0)
	    {
		/* first time seen in scan1-col */
		ep->rdeg = ep->ncolsleft + rdeg0 ;
		DEBUG6 (("e "ID" First seen: rdeg: "ID" ", e, ep->rdeg-rdeg0)) ;
		ASSERT (ep->ncolsleft > 0 && ep->nrowsleft > 0) ;
	    }

	    ep->rdeg-- ;	/* decrement external row degree */
	    DEBUG6 (("e "ID" New ext row degree: "ID"\n", e, ep->rdeg-rdeg0)) ;

	    if (ep->rdeg == rdeg0 && ep->next == EMPTY)
	    {
		/* A new LUson or Lson has been found */
		ep->next = son_list ;
		son_list = e ;
	    }

	    ASSERT (ep->rdeg >= rdeg0) ;
	    *tp2++ = *tp ;	/* leave the tuple in the list */
	}
	Col_tlen [col] = tp2 - tp1 ;
    }

    /* ---------------------------------------------------------------------- */
    /* assemble new sons via full scans */
    /* ---------------------------------------------------------------------- */

    next = EMPTY ;

    for (e = son_list ; e > 0 ; e = next)
    {
	ASSERT (e > 0 && e <= Work->nel && E [e]) ;
	p = Memory + E [e] ;
	DEBUG2 (("New son: "ID"\n", e)) ;
#ifndef NDEBUG
	UMF_dump_element (Numeric, Work, e, FALSE) ;
#endif
	GET_ELEMENT (ep, p, Cols, Rows, ncols, nrows, C) ;
	nrowsleft = ep->nrowsleft ;
	ncolsleft = ep->ncolsleft ;
	next = ep->next ;
	ep->next = EMPTY ;

	extrdeg = (ep->rdeg < rdeg0) ? ncolsleft : (ep->rdeg - rdeg0) ;
	extcdeg = (ep->cdeg < cdeg0) ? nrowsleft : (ep->cdeg - cdeg0) ;
	ncols_to_assemble = ncolsleft - extrdeg ;
	nrows_to_assemble = nrowsleft - extcdeg ;

	if (extrdeg == 0 && extcdeg == 0)
	{

	    /* -------------------------------------------------------------- */
	    /* this is an LUson - assemble an entire contribution block */
	    /* -------------------------------------------------------------- */

	    DEBUG6 (("LUson found: "ID"\n", e)) ;

	    if (nrows == nrowsleft)
	    {
		/* ---------------------------------------------------------- */
		/* no rows assembled out of this LUson yet */
		/* ---------------------------------------------------------- */

		/* compute the compressed column offset vector*/
		/* [ use Wm [0..nrows-1] for offsets */
		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    Row_degree [row] -= ncolsleft ;
		    Wm [i] = Frpos [row] ;
		}

		if (ncols == ncolsleft)
		{
		    /* ------------------------------------------------------ */
		    /* no rows or cols assembled out of LUson yet */
		    /* ------------------------------------------------------ */

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			Col_degree [col] -= nrowsleft ;
			Fcol = Fx + Fcpos [col] ;
			for (i = 0 ; i < nrows ; i++)
			{
			    /* Fcol [Wm [i]] += C [i] ; */
			    ASSEMBLE (Fcol [Wm [i]], C [i]) ;
			}
			C += nrows ;
		    }
		}
		else
		{
		    /* ------------------------------------------------------ */
		    /* only cols have been assembled out of LUson */
		    /* ------------------------------------------------------ */

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			if (col >= 0)
			{
			    Col_degree [col] -= nrowsleft ;
			    Fcol = Fx + Fcpos [col] ;
			    for (i = 0 ; i < nrows ; i++)
			    {
				/* Fcol [Wm [i]] += C [i] ; */
				ASSEMBLE (Fcol [Wm [i]], C [i]) ;
			    }
			}
			C += nrows ;
		    }
		}
		/* ] done using Wm [0..nrows-1] for offsets */
	    }
	    else
	    {
		/* ---------------------------------------------------------- */
		/* some rows have been assembled out of this LUson */
		/* ---------------------------------------------------------- */

		/* compute the compressed column offset vector*/
		/* [ use Woo,Wm [0..nrowsleft-1] for offsets */
		ngetrows = 0 ;
		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    if (row >= 0)
		    {
			Row_degree [row] -= ncolsleft ;
			Woo [ngetrows] = i ;
			Wm [ngetrows++] = Frpos [row] ;
		    }
		}
		ASSERT (ngetrows == nrowsleft) ;

		if (ncols == ncolsleft)
		{
		    /* ------------------------------------------------------ */
		    /* only rows have been assembled out of this LUson */
		    /* ------------------------------------------------------ */

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			Col_degree [col] -= nrowsleft ;
			Fcol = Fx + Fcpos [col] ;
			for (i = 0 ; i < nrowsleft ; i++)
			{
			    /* Fcol [Wm [i]] += C [Woo [i]] ; */
			    ASSEMBLE (Fcol [Wm [i]], C [Woo [i]]) ;
			}
			C += nrows ;
		    }

		}
		else
		{
		    /* ------------------------------------------------------ */
		    /* both rows and columns have been assembled out of LUson */
		    /* ------------------------------------------------------ */

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			if (col >= 0)
			{
			    Col_degree [col] -= nrowsleft ;
			    Fcol = Fx + Fcpos [col] ;
			    for (i = 0 ; i < nrowsleft ; i++)
			    {
				/* Fcol [Wm [i]] += C [Woo [i]] ; */
				ASSEMBLE (Fcol [Wm [i]], C [Woo [i]]) ;
			    }
			}
			C += nrows ;
		    }
		}
		/* ] done using Woo,Wm [0..nrowsleft-1] */
	    }

	    /* deallocate the element: remove from ordered list */
	    UMF_mem_free_tail_block (Numeric, E [e]) ;
	    E [e] = 0 ;

	}
	else if (extcdeg == 0)
	{

	    /* -------------------------------------------------------------- */
	    /* this is a Uson - assemble all possible columns */
	    /* -------------------------------------------------------------- */

	    DEBUG6 (("New USON: "ID"\n", e)) ;
	    ASSERT (extrdeg > 0) ;

	    DEBUG6 (("New uson "ID" cols to do "ID"\n", e, ncols_to_assemble)) ;

	    if (ncols_to_assemble > 0)
	    {

		if (nrows == nrowsleft)
		{
		    /* ------------------------------------------------------ */
		    /* no rows have been assembled out of this Uson yet */
		    /* ------------------------------------------------------ */

		    /* compute the compressed column offset vector */
		    /* [ use Wm [0..nrows-1] for offsets */
		    for (i = 0 ; i < nrows ; i++)
		    {
			row = Rows [i] ;
			ASSERT (row >= 0 && row < n_row) ;
			Row_degree [row] -= ncols_to_assemble ;
			Wm [i] = Frpos [row] ;
		    }
		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			if ((col >= 0) && (Fcpos [col] >= 0))
			{
			    Col_degree [col] -= nrowsleft ;
			    Fcol = Fx + Fcpos [col] ;
			    for (i = 0 ; i < nrows ; i++)
			    {
				/* Fcol [Wm [i]] += C [i] ; */
				ASSEMBLE (Fcol [Wm [i]], C [i]) ;
			    }
			    /* flag the column as assembled from Uson */
			    Cols [j] = EMPTY ;
			}
			C += nrows ;
		    }
		    /* ] done using Wm [0..nrows-1] for offsets */
		}
		else
		{
		    /* ------------------------------------------------------ */
		    /* some rows have been assembled out of this Uson */
		    /* ------------------------------------------------------ */

		    /* compute the compressed column offset vector*/
		    /* [ use Woo,Wm [0..nrows-1] for offsets */
		    ngetrows = 0 ;
		    for (i = 0 ; i < nrows ; i++)
		    {
			row = Rows [i] ;
			if (row >= 0)
			{
			    Row_degree [row] -= ncols_to_assemble ;
			    ASSERT (row < n_row && Frpos [row] >= 0) ;
			    Woo [ngetrows] = i ;
			    Wm [ngetrows++] = Frpos [row] ;
			}
		    }
		    ASSERT (ngetrows == nrowsleft) ;

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			if ((col >= 0) && (Fcpos [col] >= 0))
			{
			    Col_degree [col] -= nrowsleft ;
			    Fcol = Fx + Fcpos [col] ;
			    for (i = 0 ; i < nrowsleft ; i++)
			    {
				/* Fcol [Wm [i]] += C [Woo [i]] ; */
				ASSEMBLE (Fcol [Wm [i]], C [Woo [i]]) ;
			    }
			    /* flag the column as assembled from Uson */
			    Cols [j] = EMPTY ;
			}
			C += nrows ;
		    }
		    /* ] done using Woo,Wm */
		}
		ep->ncolsleft = extrdeg ;
	    }

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* this is an Lson - assemble all possible rows */
	    /* -------------------------------------------------------------- */

	    DEBUG6 (("New LSON: "ID"\n", e)) ;
	    ASSERT (extrdeg == 0 && extcdeg > 0) ;

	    DEBUG6 (("New lson "ID" rows to do "ID"\n", e, nrows_to_assemble)) ;

	    if (nrows_to_assemble > 0)
	    {

		/* compute the compressed column offset vector */
		/* [ use Woo,Wm [0..nrows-1] for offsets */
		ngetrows = 0 ;
		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    if ((row >= 0) && (Frpos [row] >= 0))
		    {
			ASSERT (row < n_row) ;
			Row_degree [row] -= ncolsleft ;
			Woo [ngetrows] = i ;
			Wm [ngetrows++] = Frpos [row] ;
			/* flag the row as assembled from the Lson */
			Rows [i] = EMPTY ;
		    }
		}
		ASSERT (nrowsleft - ngetrows == extcdeg) ;
		ASSERT (ngetrows == nrows_to_assemble) ;

		if (ncols == ncolsleft)
		{
		    /* ------------------------------------------------------ */
		    /* no columns assembled out this Lson yet */
		    /* ------------------------------------------------------ */

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			ASSERT (col >= 0 && col < n_col) ;
			Col_degree [col] -= nrows_to_assemble ;
			Fcol = Fx + Fcpos [col] ;
			for (i = 0 ; i < nrows_to_assemble ; i++)
			{
			    /* Fcol [Wm [i]] += C [Woo [i]] ; */
			    ASSEMBLE (Fcol [Wm [i]], C [Woo [i]]) ;
			}
			C += nrows ;
		    }
		}
		else
		{
		    /* ------------------------------------------------------ */
		    /* some columns have been assembled out of this Lson */
		    /* ------------------------------------------------------ */

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			ASSERT (col < n_col) ;
			if (col >= 0)
			{
			    Col_degree [col] -= nrows_to_assemble ;
			    Fcol = Fx + Fcpos [col] ;
			    for (i = 0 ; i < nrows_to_assemble ; i++)
			    {
				/* Fcol [Wm [i]] += C [Woo [i]] ; */
				ASSEMBLE (Fcol [Wm [i]], C [Woo [i]]) ;
			    }
			}
			C += nrows ;
		    }
		}

		/* ] done using Woo,Wm */

		ep->nrowsleft = extcdeg ;
	    }
	}
    }

    /* Note that garbage collection, and build tuples */
    /* both destroy the son list. */

    /* ] son_list now empty */

    /* ---------------------------------------------------------------------- */
    /* If frontal matrix extended, assemble old L/Usons from new rows/cols */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* SCAN2-row:  assemble rows of old Lsons from the new rows */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    DEBUG7 (("Current frontal matrix: (prior to scan2-row)\n")) ;
    UMF_dump_current_front (Numeric, Work, TRUE) ;
#endif

    if (Work->do_scan2row)
    {

	scan_pivrow = Work->scan_pivrow ;

	for (i2 = Work->fscan_row ; i2 < fnrows || scan_pivrow ; )
	{
	    /* Get a row */
	    if (scan_pivrow)
	    {
		/* pivot row is new to this front.  Scan it */
		row = Work->pivrow ;
		scan_pivrow = FALSE ;
	    }
	    else
	    {
		row = Frows [i2++] ;
	    }
	    ASSERT (row >= 0 && row < n_row) ;

#ifndef NDEBUG
	    DEBUG6 (("SCAN2-row: "ID"\n", row)) ;
	    UMF_dump_rowcol (0, Numeric, Work, row, FALSE) ;
#endif

	    ASSERT (NON_PIVOTAL_ROW (row)) ;
	    tpi = Row_tuples [row] ;
	    if (!tpi) continue ;
	    tp = (Tuple *) (Memory + tpi) ;
	    tp1 = tp ;
	    tp2 = tp ;
	    tpend = tp + Row_tlen [row] ;
	    for ( ; tp < tpend ; tp++)
	    {
		e = tp->e ;
		ASSERT (e > 0 && e <= Work->nel) ;
		if (!E [e]) continue ;	/* element already deallocated */
		f = tp->f ;
		p = Memory + E [e] ;
		ep = (Element *) p ;
		p += UNITS (Element, 1) ;
		Cols = (Int *) p ;
		Rows = Cols + ep->ncols ;
		if (Rows [f] == EMPTY) continue ;    /* row already assembled */
		ASSERT (row == Rows [f] && row >= 0 && row < n_row) ;

		if (ep->rdeg == rdeg0)
		{
		    /* ------------------------------------------------------ */
		    /* this is an old Lson - assemble just one row */
		    /* ------------------------------------------------------ */

		    /* flag the row as assembled from the Lson */
		    Rows [f] = EMPTY ;

		    nrows = ep->nrows ;
		    ncols = ep->ncols ;
		    p += UNITS (Int, ncols + nrows) ;
		    C = ((Entry *) p) + f ;

		    DEBUG6 (("Old LSON: "ID"\n", e)) ;
#ifndef NDEBUG
		    UMF_dump_element (Numeric, Work, e, FALSE) ;
#endif

		    nrowsleft = ep->nrowsleft ;
		    ncolsleft = ep->ncolsleft ;

		    Frow = Fx + Frpos [row] ;
		    DEBUG6 (("LSON found (in scan2-row): "ID"\n", e)) ;

		    Row_degree [row] -= ncolsleft ;

		    if (ncols == ncolsleft)
		    {
			/* -------------------------------------------------- */
			/* no columns assembled out this Lson yet */
			/* -------------------------------------------------- */

			for (j = 0 ; j < ncols ; j++)
			{
			    col = Cols [j] ;
			    ASSERT (col >= 0 && col < n_col) ;
			    Col_degree [col] -- ;
			    /* Frow [Fcpos [col]] += *C ; */
			    ASSEMBLE (Frow [Fcpos [col]], *C) ;
			    C += nrows ;
			}
		    }
		    else
		    {
			/* -------------------------------------------------- */
			/* some columns have been assembled out of this Lson */
			/* -------------------------------------------------- */

			for (j = 0 ; j < ncols ; j++)
			{
			    col = Cols [j] ;
			    if (col >= 0)
			    {
				ASSERT (col < n_col) ;
				Col_degree [col] -- ;
				/* Frow [Fcpos [col]] += *C ; */
				ASSEMBLE (Frow [Fcpos [col]], *C) ;
			    }
			    C += nrows ;
			}
		    }
		    ep->nrowsleft-- ;
		    ASSERT (ep->nrowsleft > 0) ;
		}
		else
		{
		    *tp2++ = *tp ;	/* leave the tuple in the list */
		}
	    }
	    Row_tlen [row] = tp2 - tp1 ;

#ifndef NDEBUG
	    DEBUG7 (("row assembled in scan2-row: "ID"\n", row)) ;
	    if (row != Work->pivrow)
	    {
		UMF_dump_rowcol (0, Numeric, Work, row, FALSE) ;
	    }
	    DEBUG7 (("Current frontal matrix: (scan 1b)\n")) ;
	    UMF_dump_current_front (Numeric, Work, TRUE) ;
#endif
	}
    }

    /* ---------------------------------------------------------------------- */
    /* SCAN2-col:  assemble columns of old Usons from the new columns */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    DEBUG7 (("Current frontal matrix: (prior to scan2-col)\n")) ;
    UMF_dump_current_front (Numeric, Work, TRUE) ;
#endif

    if (Work->do_scan2col)
    {

	scan_pivcol = Work->scan_pivcol ;

	for (j2 = Work->fscan_col ; j2 < fncols || scan_pivcol ; )
	{
	    /* Get a column */
	    if (scan_pivcol)
	    {
		/* pivot col is new to this front.  Scan it */
		col = Work->pivcol ;
		scan_pivcol = FALSE ;
	    }
	    else
	    {
		col = Fcols [j2++] ;
	    }
	    ASSERT (col >= 0 && col < n_col) ;

	    DEBUG6 (("SCAN2-col: "ID"\n", col)) ;
#ifndef NDEBUG
	    UMF_dump_rowcol (1, Numeric, Work, col, FALSE) ;
#endif

	    ASSERT (NON_PIVOTAL_COL (col)) ;
	    tpi = Col_tuples [col] ;
	    if (!tpi) continue ;
	    tp = (Tuple *) (Memory + tpi) ;
	    tp1 = tp ;
	    tp2 = tp ;
	    tpend = tp + Col_tlen [col] ;
	    for ( ; tp < tpend ; tp++)
	    {

		e = tp->e ;
		ASSERT (e > 0 && e <= Work->nel) ;
		if (!E [e]) continue ;	/* element already deallocated */
		f = tp->f ;
		p = Memory + E [e] ;
		ep = (Element *) p ;
		p += UNITS (Element, 1) ;
		Cols = (Int *) p ;
		if (Cols [f] == EMPTY)
		{
		    continue ;	/* column already assembled */
		}
		ASSERT (col == Cols [f] && col >= 0 && col < n_col) ;

		if (ep->cdeg == cdeg0)
		{

		    /* ------------------------------------------------------ */
		    /* this is an old Uson - assemble just one column */
		    /* ------------------------------------------------------ */

		    /* flag the column as assembled from the Uson */
		    Cols [f] = EMPTY ;

		    nrows = ep->nrows ;
		    ncols = ep->ncols ;
		    Rows = Cols + ncols ;
		    p += UNITS (Int, ep->ncols + nrows) ;
		    C = ((Entry *) p) + f * nrows ;

		    DEBUG6 (("Old USON: "ID"\n", e)) ;
#ifndef NDEBUG
		    UMF_dump_element (Numeric, Work, e, FALSE) ;
#endif

		    nrowsleft = ep->nrowsleft ;
		    ncolsleft = ep->ncolsleft ;

		    Fcol = Fx + Fcpos [col] ;
		    DEBUG6 (("USON found (in scan2-col): "ID"\n", e)) ;

		    Col_degree [col] -= nrowsleft ;

		    if (nrows == nrowsleft)
		    {
			/* -------------------------------------------------- */
			/* no rows assembled out of this Uson yet */
			/* -------------------------------------------------- */

			for (i = 0 ; i < nrows ; i++)
			{
			    row = Rows [i] ;
			    ASSERT (row >= 0 && row < n_row) ;
			    Row_degree [row]-- ;
			    /* Fcol [Frpos [row]] += C [i] ; */
			    ASSEMBLE (Fcol [Frpos [row]], C [i]) ;
			}
		    }
		    else
		    {
			/* -------------------------------------------------- */
			/* some rows have been assembled out of this Uson */
			/* -------------------------------------------------- */

			for (i = 0 ; i < nrows ; i++)
			{
			    row = Rows [i] ;
			    if (row >= 0)
			    {
				ASSERT (row < n_row) ;
				Row_degree [row]-- ;
				/* Fcol [Frpos [row]] += C [i] ; */
				ASSEMBLE (Fcol [Frpos [row]], C [i]) ;
			    }
			}
		    }
		    ep->ncolsleft-- ;
		    ASSERT (ep->ncolsleft > 0) ;
		}
		else
		{
		    *tp2++ = *tp ;	/* leave the tuple in the list */
		}
	    }
	    Col_tlen [col] = tp2 - tp1 ;

#ifndef NDEBUG
	    DEBUG7 (("Column assembled in scan2-col: "ID"\n", col)) ;
	    if (col != Work->pivcol)
	    {
		UMF_dump_rowcol (1, Numeric, Work, col, FALSE) ;
	    }
	    DEBUG7 (("Current frontal matrix: after scan2-col\n")) ;
	    UMF_dump_current_front (Numeric, Work, TRUE) ;
#endif

	}

    }

    /* ---------------------------------------------------------------------- */
    /* done.  the remainder of this routine is used only when in debug mode */
    /* ---------------------------------------------------------------------- */



#ifndef NDEBUG

    /* ---------------------------------------------------------------------- */
    /* when debugging: make sure the assembly did everything that it could */
    /* ---------------------------------------------------------------------- */

    DEBUG3 (("::Assemble done\n")) ;

    scan_pivrow = TRUE ;

    for (i2 = 0 ; i2 < fnrows || scan_pivrow ; )
    {
	/* Get a row */
	if (scan_pivrow)
	{
	    /* pivot row is new to this front.  Scan it */
	    row = Work->pivrow ;
	    scan_pivrow = FALSE ;
	}
	else
	{
	    row = Frows [i2++] ;
	}
	ASSERT (row >= 0 && row < n_row) ;

	DEBUG6 (("DEBUG SCAN 1: "ID"\n", row)) ;
	UMF_dump_rowcol (0, Numeric, Work, row, TRUE) ;

	ASSERT (NON_PIVOTAL_ROW (row)) ;
	tpi = Row_tuples [row] ;
	if (!tpi) continue ;
	tp = (Tuple *) (Memory + tpi) ;
	tpend = tp + Row_tlen [row] ;
	for ( ; tp < tpend ; tp++)
	{
	    e = tp->e ;
	    ASSERT (e > 0 && e <= Work->nel) ;
	    if (!E [e]) continue ;	/* element already deallocated */
	    f = tp->f ;
	    p = Memory + E [e] ;
	    ep = (Element *) p ;
	    p += UNITS (Element, 1) ;
	    Rows = ((Int *) p) + ep->ncols ;
	    if (Rows [f] == EMPTY) continue ;	/* row already assembled */
	    ASSERT (row == Rows [f]) ;
	    extrdeg = (ep->rdeg < rdeg0) ? ep->ncolsleft : (ep->rdeg - rdeg0) ;
	    extcdeg = (ep->cdeg < cdeg0) ? ep->nrowsleft : (ep->cdeg - cdeg0) ;
	    DEBUG6 ((
		"e "ID" After assembly ext row deg: "ID" ext col degree "ID"\n",
		e, extrdeg, extcdeg)) ;

	    /* no Lsons in any row */
	    ASSERT (extrdeg > 0) ;

	    /* Uson external row degree is = number of cols left */
	    ASSERT (IMPLIES (extcdeg == 0, extrdeg == ep->ncolsleft)) ;
	}
    }

    /* ---------------------------------------------------------------------- */

    scan_pivcol = TRUE ;

    for (j2 = 0 ; j2 < fncols || scan_pivcol ; )
    {
	/* Get a column */
	if (scan_pivcol)
	{
	    /* pivot col is new to this front.  Scan it */
	    col = Work->pivcol ;
	    scan_pivcol = FALSE ;
	}
	else
	{
	    col = Fcols [j2++] ;
	}
	ASSERT (col >= 0 && col < n_col) ;

	DEBUG6 (("DEBUG SCAN 2: "ID"\n", col)) ;
	UMF_dump_rowcol (1, Numeric, Work, col, TRUE) ;

	ASSERT (NON_PIVOTAL_COL (col)) ;
	tpi = Col_tuples [col] ;
	if (!tpi) continue ;
	tp = (Tuple *) (Memory + tpi) ;
	tpend = tp + Col_tlen [col] ;
	for ( ; tp < tpend ; tp++)
	{
	    e = tp->e ;
	    ASSERT (e > 0 && e <= Work->nel) ;
	    if (!E [e]) continue ;	/* element already deallocated */
	    f = tp->f ;
	    p = Memory + E [e] ;
	    ep = (Element *) p ;
	    p += UNITS (Element, 1) ;
	    Cols = (Int *) p ;
	    if (Cols [f] == EMPTY) continue ;	/* column already assembled */
	    ASSERT (col == Cols [f]) ;
	    extrdeg = (ep->rdeg < rdeg0) ? ep->ncolsleft : (ep->rdeg - rdeg0) ;
	    extcdeg = (ep->cdeg < cdeg0) ? ep->nrowsleft : (ep->cdeg - cdeg0) ;
	    DEBUG6 (("e "ID" After assembly ext col deg: "ID"\n", e, extcdeg)) ;

	    /* no Usons in any column */
	    ASSERT (extcdeg > 0) ;

	    /* Lson external column degree is = number of rows left */
	    ASSERT (IMPLIES (extrdeg == 0, extcdeg == ep->nrowsleft)) ;
	}
    }

#endif /* NDEBUG */

}

