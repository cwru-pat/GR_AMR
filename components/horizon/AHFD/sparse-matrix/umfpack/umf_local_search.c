/* ========================================================================== */
/* === UMF_local_search ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Perform pivot search to find pivot row and pivot column.

    The pivot column is selected from the candidate set (a supercolumn from
    colamd), and then removed from that set.

    Called by kernel.

    Modifies front, but keeps current Frows and Fcols unchanged (new entries
    are appended).  This can be undone if the current front must be written
    out after the pivot is found.  The numerical values of the front are
    unchanged.

    The current frontal matrix might be empty (Fx not allocated, Work->fnrows
    and Work->fncols are both zero).

    Returns UMFPACK_OK if successful, or UMFPACK_WARNING_singular_matrix, or
    UMFPACK_ERROR_different_pattern if not.

*/

/* ========================================================================== */

#include "umf_internal.h"
#include "umf_row_search.h"
#include "umf_mem_free_tail_block.h"

#define IN 0
#define OUT 1

#define IN_IN 0
#define IN_OUT 1
#define OUT_IN 2
#define OUT_OUT 3

/* ========================================================================== */

GLOBAL Int UMF_local_search
(
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Entry *Fx, *Fl, *Fu, *Fs, *C, *Fd ;
    double relax, relax2, relax3 ;
    Int pos, nrows, *Cols, *Rows, e, f, jmax, status, max_cdeg, fnzeros,
	j, col, i, row, cdeg [2], rdeg [2][2], fnpiv, nothing [2], new_LUsize,
	pivrow [2][2], pivcol [2], *Wp, *Fcpos, *Frpos, new_fnzeros,
	*Wm, *Wio, *Woi, *Woo, *Frows, *Fcols, fnrows, fncols, *E,
	deg, nr_in, nc, src, dest, thiscost, bestcost, nr_out, *Wcol, do_update,
	extra_cols, extra_rows, extra_zeros, relaxed_front, do_extend,
	fnrows_max, fncols_max, tpi, *Col_tuples, *Col_degree, *Col_tlen,
	jj, jcand [2], freebie [2], did_rowmerge ;
    Unit *Memory, *p ;
    Tuple *tp, *tpend, *tp1, *tp2 ;
    Element *ep ;

#ifndef NDEBUG
    Int debug_ok, n_row, n_col, *Row_degree ;
    Row_degree = Numeric->Rperm ;	/* for NON_PIVOTAL_ROW macro only */
#endif

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    Memory = Numeric->Memory ;
    E = Work->E ;
    Col_degree = Numeric->Cperm ;

    Col_tuples = Numeric->Lip ;
    Col_tlen   = Numeric->Lilen ;

    Wp = Work->Wp ;
    Wm = Work->Wm ;
    Woi = Work->Woi ;
    Wio = Work->Wio ;
    Woo = Work->Woo ;	/* size at least fncols_max */
    Fx = Work->Fx ;
    Fcpos = Work->Fcpos ;
    Frpos = Work->Frpos ;
    Frows = Work->Frows ;
    Fcols = Work->Fcols ;
    fnrows = Work->fnrows ;
    fncols = Work->fncols ;
    fnrows_max = Work->fnrows_max ;
    fncols_max = Work->fncols_max ;
    fnpiv = Work->fnpiv ;
    nothing [0] = EMPTY ;
    nothing [1] = EMPTY ;
    relax = Numeric->relax ;
    relax2 = Numeric->relax2 ;
    relax3 = Numeric->relax3 ;
    fnzeros = Work->fnzeros ;
    new_fnzeros = fnzeros ;
    jj = EMPTY ;

    /* The pivot column degree cannot exceed max_cdeg */
    max_cdeg = fnrows_max - fnpiv ;

#ifndef NDEBUG
    n_row = Work->n_row ;
    n_col = Work->n_col ;
    DEBUG2 (("LOCAL SEARCH:  current frontal matrix: # candidates: "ID"\n",
	Work->ncand)) ;
    UMF_dump_current_front (Numeric, Work, TRUE) ;
    if (UMF_debug > 0 || MAX (n_row, n_col) < 1000)
    {
	for (i = 0 ; i < MAX (n_row, n_col) ; i++)
	{
	    ASSERT (Wp [i] < 0) ;
	    ASSERT (Wp [i] > Work->Wpflag) ;
	}
    }
    ASSERT (Work->Wpflag < EMPTY) ;
    DEBUG2 (("candidates: ")) ;
    for (j = 0 ; j < Work->ncand ; j++) DEBUG2 ((ID" ", Work->Candidates [j]));
    DEBUG2 (("\n")) ;
#endif

    /* ====================================================================== */
    /* === PIVOT SEARCH ===================================================== */
    /* ====================================================================== */

    /* initialize */

    pivcol [IN] = EMPTY ;
    pivcol [OUT] = EMPTY ;

    /* Int_MAX is defined in umfpack.h */
    cdeg [IN] = Int_MAX ;
    cdeg [OUT] = Int_MAX ;

    pivrow [IN][IN] = EMPTY ;
    pivrow [IN][OUT] = EMPTY ;
    pivrow [OUT][IN] = EMPTY ;
    pivrow [OUT][OUT] = EMPTY ;

    rdeg [IN][IN] = Int_MAX ;
    rdeg [IN][OUT] = Int_MAX ;
    rdeg [OUT][IN] = Int_MAX ;
    rdeg [OUT][OUT] = Int_MAX ;

    freebie [IN] = FALSE ;
    freebie [OUT] = FALSE ;

    Work->pivot_case = EMPTY ;
    bestcost = EMPTY ;

    nr_out = EMPTY ;
    nr_in = EMPTY ;

    jcand [IN] = EMPTY ;
    jcand [OUT] = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* find shortest column in the front, and shortest column not in the */
    /* front, from the candidate pivot column set */
    /* ---------------------------------------------------------------------- */

    /* If there are too many candidates, then only look at the first */
    /* MAX_CANDIDATES of them.   Otherwise, if there are O(n) candidates, */
    /* this code could take O(n^2) time. */

    jmax = MIN (MAX_CANDIDATES, Work->ncand) ;

#ifndef NDEBUG
    DEBUG2 (("Max candidates "ID", Work->ncand "ID" jmax "ID"\n",
	MAX_CANDIDATES, Work->ncand, jmax)) ;
#endif

    /* get the first candidate */

    col = Work->Candidates [0] ;
    DEBUG3 (("Pivot column candidate: "ID" j = "ID"\n", col, j)) ;
    ASSERT (col >= 0 && col < n_col) ;
    deg = Col_degree [col] ;

#ifndef NDEBUG
    DEBUG3 (("Pivot column candidate: "ID" cost: "ID"  Fcpos[col] "ID"\n",
	col, deg, Fcpos [col])) ;
    UMF_dump_rowcol (1, Numeric, Work, col, TRUE) ;
#endif

    if (Fcpos [col] >= 0)
    {
	/* best column in front, so far */
	pivcol [IN] = col ;
	cdeg [IN] = deg ;
	jcand [IN] = 0 ;
    }
    else
    {
	/* best column not in front, so far */
	pivcol [OUT] = col ;
	cdeg [OUT] = deg ;
	jcand [OUT] = 0 ;
    }

    /* look at the rest of the candidates */

    for (j = 1 ; j < jmax ; j++)
    {
	col = Work->Candidates [j] ;
	DEBUG3 (("Pivot column candidate: "ID" j = "ID"\n", col, j)) ;
	ASSERT (col >= 0 && col < n_col) ;
	deg = Col_degree [col] ;
#ifndef NDEBUG
	DEBUG3 (("Pivot column candidate: "ID" cost: "ID"  Fcpos[col] "ID"\n",
		col, deg, Fcpos [col])) ;
	UMF_dump_rowcol (1, Numeric, Work, col, TRUE) ;
#endif
	if (Fcpos [col] >= 0)
	{
#ifndef NDEBUG
	    Int fs ;
	    fs = Fcpos [col] / fnrows_max ;
	    ASSERT (fs >= 0 && fs < fncols) ;
#endif
	    if (deg < cdeg [IN] || (deg == cdeg [IN] && col < pivcol [IN]))
	    {
		/* best column in front, so far */
		pivcol [IN] = col ;
		cdeg [IN] = deg ;
		jcand [IN] = j ;
	    }
	}
	else
	{
	    if (deg < cdeg [OUT] || (deg == cdeg [OUT] && col < pivcol [OUT]))
	    {
		/* best column not in front, so far */
		pivcol [OUT] = col ;
		cdeg [OUT] = deg ;
		jcand [OUT] = j ;
	    }
	}
    }

    DEBUG2 (("Pivcol in "ID" out "ID"\n", pivcol [IN], pivcol [OUT])) ;
    ASSERT ((pivcol [IN] >= 0 && pivcol [IN] < n_col)
	|| (pivcol [OUT] >= 0 && pivcol [OUT] < n_col)) ;

    cdeg [IN] = EMPTY ;
    cdeg [OUT] = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* construct candidate column in front, and search for pivot rows */
    /* ---------------------------------------------------------------------- */

    if (pivcol [IN] != EMPTY)
    {

#ifndef NDEBUG
	DEBUG2 (("col[IN] column "ID" in front at position = "ID"\n",
		pivcol [IN], Fcpos [pivcol [IN]])) ;
	UMF_dump_rowcol (1, Numeric, Work, pivcol [IN], TRUE) ;
#endif

	/* the only way we can have a pivcol[IN] is if the front is not empty */
	ASSERT (fnrows > 0 && fncols > 0) ;

	Fs = Fx + Fcpos [pivcol [IN]] ;

	/* ------------------------------------------------------------------ */
	/* update column in front (permanent) */
	/* ------------------------------------------------------------------ */

	DEBUG6 (("Update pivot column:\n")) ;
	Fl = Fx + (fncols_max - fnpiv) * fnrows_max ;
	Fu = Fs + (fnrows_max - fnpiv) ;

#ifdef USE_NO_BLAS

	/* no BLAS available - use plain C code instead */
	for (j = 0 ; j < fnpiv ; j++)
	{
	    Entry Fuj ;
	    Fuj = Fu [j] ;
	    if (IS_NONZERO (Fuj))
	    {
		for (i = 0 ; i < fnrows ; i++)
		{
		    /* Fs [i] -= Fuj * Fl [i+j*fnrows_max] ; */
		    MULT_SUB (Fs [i], Fuj, Fl [i+j*fnrows_max]) ;
		}
	    }
	}

#else

	BLAS_GEMV_COL (fnrows, fnpiv, Fl, Fu, Fs, fnrows_max) ;

#endif /* USE_NO_BLAS */

	/* zero out the column of U, in case this doesn't become pivot column */
	for (j = 0 ; j < fnpiv ; j++)
	{
	    CLEAR (Fu [j]) ;
	}

	/* ------------------------------------------------------------------ */

#ifndef NDEBUG
	DEBUG6 (("Fs after update: fnrows="ID"\n", fnrows)) ;
	DEBUG6 ((" Work->fnpiv="ID" \n", fnpiv)) ;
	for (i = 0 ; i < fnrows ; i++)
	{
	    DEBUG6 ((ID" "ID" "ID, i, Frows [i], Frpos [Frows [i]])) ;
	    EDEBUG6 (Fs [i]) ;
	    DEBUG6 (("\n")) ;
	}
#endif

	/* ------------------------------------------------------------------ */
	/* construct the candidate column in the front */
	/* ------------------------------------------------------------------ */

	cdeg [IN] = fnrows ;

#ifndef NDEBUG
	/* check Frpos */
	DEBUG5 (("COL ASSEMBLE: cdeg "ID"\nREDUCE COL in "ID" max_cdeg "ID"\n",
		cdeg [IN], pivcol [IN], max_cdeg)) ;
	for (i = 0 ; i < cdeg [IN] ; i++)
	{
	    row = Frows [i] ;
	    ASSERT (row >= 0 && row < n_row) ;
	    ASSERT (Frpos [row] == i) ;
	}
	if (UMF_debug > 0 || n_row < 1000)
	{
	    Int cnt = cdeg [IN] ;
	    for (row = 0 ; row < n_row ; row++)
	    {
		if (Frpos [row] == EMPTY) cnt++ ;
	    }
	    ASSERT (cnt == n_row) ;
	}
#endif

	ASSERT (pivcol [IN] >= 0 && pivcol [IN] < n_col) ;
	ASSERT (NON_PIVOTAL_COL (pivcol [IN])) ;

	tpi = Col_tuples [pivcol [IN]] ;
	if (tpi)
	{
	    tp = (Tuple *) (Memory + tpi) ;
	    tp1 = tp ;
	    tp2 = tp ;
	    tpend = tp + Col_tlen [pivcol [IN]] ;
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
		if (Cols [f] == EMPTY) continue ; /* column already assembled */
		ASSERT (pivcol [IN] == Cols [f]) ;

		Rows = Cols + ep->ncols ;
		nrows = ep->nrows ;
		p += UNITS (Int, ep->ncols + nrows) ;
		C = ((Entry *) p) + f * nrows ;

		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    if (row >= 0) /* skip this if already gone from element */
		    {
			ASSERT (row < n_row) ;
			pos = Frpos [row] ;
			if (pos < 0)
			{
			    if (cdeg [IN] >= max_cdeg)
			    {
				return (UMFPACK_ERROR_different_pattern) ;
			    }
			    /* new entry in the pattern - save Frpos */
			    Frpos [row] = cdeg [IN] ;
			    Frows [cdeg [IN]] = row ;

			    /* this will be discarded */
			    Fs [cdeg [IN]++] = C [i] ;
			}
			else
			{
			    /* entry already in pattern - sum the values */
			    /* Fs [pos] += C [i] ; */
			    ASSEMBLE (Fs [pos], C [i]) ;
			    if (pos < fnrows)
			    {
				/* just did a permanent assembly of a single */
				/* entry into the current front */
				CLEAR (C [i]) ;
			    }
			}
		    }
		}

		*tp2++ = *tp ;	/* leave the tuple in the list */

	    }
	    Col_tlen [pivcol [IN]] = tp2 - tp1 ;
	}

	/* ------------------------------------------------------------------ */

#ifndef NDEBUG
	DEBUG4 (("Reduced column: cdeg in  "ID"\n", cdeg [IN])) ;
	for (i = 0 ; i < cdeg [IN] ; i++)
	{
	    DEBUG6 ((" "ID" "ID" "ID, i, Frows [i], Frpos [Frows [i]])) ;
	    EDEBUG6 (Fs [i]) ;
	    DEBUG6 (("\n")) ;
	    ASSERT (i == Frpos [Frows [i]]) ;
	}
	ASSERT (cdeg [IN] + fnpiv <= fnrows_max) ;
#endif

	/* ------------------------------------------------------------------ */
	/* cdeg [IN] is now the exact degree of this column */
	/* ------------------------------------------------------------------ */

	nr_in = cdeg [IN] - fnrows ;

	/* since there are no 0-by-x fronts, if there is a pivcol [IN] the */
	/* front must have at least one row. */
	ASSERT (cdeg [IN] > 0) ;

	/* new degree of pivcol [IN], excluding current front is nr_in */
	/* column expands by nr_in rows */

	/* ------------------------------------------------------------------ */
	/* search for two candidate pivot rows */
	/* ------------------------------------------------------------------ */

	/* for the IN_IN pivot row (if any), */
	/* extend the pattern in place, using Fcols */
	status = UMF_row_search (Numeric, Work, Symbolic, cdeg [IN], Frows,
	    pivrow [IN], rdeg [IN], Fcols, Wio, nothing, Fs, pivcol [IN],
	    freebie) ;
	ASSERT (!freebie [IN] && !freebie [OUT]) ;

	/* ------------------------------------------------------------------ */
	/* fatal error if matrix pattern has changed since symbolic analysis */
	/* ------------------------------------------------------------------ */

	if (status == UMFPACK_ERROR_different_pattern)
	{
	    return (UMFPACK_ERROR_different_pattern) ;
	}

	/* ------------------------------------------------------------------ */
	/* we now must have a structural pivot */
	/* ------------------------------------------------------------------ */

	/* Since the pivcol[IN] exists, there must be at least one row in the */
	/* current frontal matrix, and so we must have found a structural */
	/* pivot.  The numerical value might be zero, of course. */

	ASSERT (status != UMFPACK_WARNING_singular_matrix) ;

	/* ------------------------------------------------------------------ */
	/* evaluate IN_IN option */
	/* ------------------------------------------------------------------ */

	if (pivrow [IN][IN] != EMPTY)
	{
	    /* the current front would become an (implicit) LUson */
	    /* cost is how much the current front would expand */

	    /* pivrow[IN][IN] candidates are not found via row merge search */

	    ASSERT (cdeg [IN] > 0) ;
	    nc = rdeg [IN][IN] - fncols ;

	    thiscost =
		/* each column in front (except pivot column) grows by nr_in: */
		(nr_in * (fncols - 1)) +
		/* new columns not in old front: */
		(nc * (cdeg [IN] - 1)) ;

	    /* no extra cost to relaxed amalgamation */

	    ASSERT (fnrows + nr_in == cdeg [IN]) ;
	    ASSERT (fncols + nc == rdeg [IN][IN]) ;

	    /* relaxed_front = ((fncols-1) + nc) * ((fnrows-1) + nr_in) ; */
	    do_extend = TRUE ;

	    DEBUG2 (("Evaluating option IN-IN:\n")) ;
	    DEBUG2 (("Work->fnzeros "ID" fnpiv "ID" nr_in "ID" nc "ID"\n",
	    	Work->fnzeros, fnpiv, nr_in, nc)) ;
	    DEBUG2 (("fncols "ID" fnrows "ID"\n", fncols, fnrows)) ;

	    /* determine if BLAS-3 update should be applied before extending. */
	    /* update if too many zero entries accumulate in the LU block */
	    fnzeros = Work->fnzeros + fnpiv * (nr_in + nc) ;

	    DEBUG2 (("fnzeros "ID"\n", fnzeros)) ;

	    new_LUsize = (fnpiv+1) * (fnrows + nr_in + fncols + nc) ;

	    DEBUG2 (("new_LUsize "ID"\n", new_LUsize)) ;
	    DEBUG2 (("relax2 %g\n", relax2)) ;

	    /* relax2 parameter uses a double relop, but ignore NaN case: */
	    do_update = (((double) fnzeros) / ((double) new_LUsize)) > relax2 ;

	    DEBUG2 (("do_update "ID"\n", do_update))

	    DEBUG2 (("option IN  IN : nr "ID" nc "ID" cost "ID"(0) relax "ID
		"\n", nr_in, nc, thiscost, do_extend)) ;

	    /* this is the best option seen so far */
	    Work->pivot_case = IN_IN ;
	    bestcost = thiscost ;

	    /* do the amalgamation and extend the front */
	    Work->do_extend = do_extend ;
	    Work->do_update = do_update ;
	    new_fnzeros = fnzeros ;

	}

	/* ------------------------------------------------------------------ */
	/* evaluate IN_OUT option */
	/* ------------------------------------------------------------------ */

	if (pivrow [IN][OUT] != EMPTY)
	{
	    /* the current front would become a Uson of the new front */

	    ASSERT (cdeg [IN] > 0) ;

	    /* must be at least one row outside the front */
	    /* (the pivrow [IN][OUT] itself) */
	    ASSERT (nr_in >= 1) ;

	    /* count columns not in current front */
	    nc = 0 ;
#ifndef NDEBUG
	    debug_ok = FALSE ;
#endif
	    for (i = 0 ; i < rdeg [IN][OUT] ; i++)
	    {
		col = Wio [i] ;
		DEBUG6 (("counting col "ID" Fcpos[] = "ID"\n", col,
		    Fcpos [col])) ;
		ASSERT (col >= 0 && col < n_col) ;
		ASSERT (NON_PIVOTAL_COL (col)) ;
		if (Fcpos [col] < 0) nc++ ;
#ifndef NDEBUG
		/* we must see the pivot column somewhere */
		if (col == pivcol [IN])
		{
		    ASSERT (Fcpos [col] >= 0) ;
		    debug_ok = TRUE ;
		}
#endif
	    }
	    ASSERT (debug_ok) ;

	    thiscost =
		/* each row in front grows by nc: */
		(nc * fnrows) +
		/* new rows not affected by front: */
		((nr_in-1) * (rdeg [IN][OUT]-1)) ;

	    /* check the cost of relaxed IN_OUT amalgamation */

	    extra_cols = ((fncols-1) + nc ) - (rdeg [IN][OUT] - 1) ;
	    ASSERT (extra_cols >= 0) ;
	    ASSERT (fncols + nc == extra_cols + rdeg [IN][OUT]) ;
	    extra_zeros = (nr_in-1) * extra_cols ;	/* symbolic fill-in */

	    ASSERT (fnrows + nr_in == cdeg [IN]) ;
	    ASSERT (fncols + nc == rdeg [IN][OUT] + extra_cols) ;

	    /* size of relaxed front: */
	    relaxed_front = ((fncols-1) + nc) * (fnrows + (nr_in-1)) ;

	    /* do relaxed amalgamation if the extra zeros are no more */
	    /* than a fraction (default 0.25) of the relaxed front */
	    /* if relax = 0: no extra zeros allowed */
	    /* if relax = +inf: always amalgamate */

	    /* relax parameter uses a double relop, but ignore NaN case: */
	    do_extend = ((double) extra_zeros)
		< (relax * (double) relaxed_front) ;

	    if (do_extend)
	    {
		/* count the cost of relaxed amalgamation */
		thiscost += extra_zeros ;

	        DEBUG2 (("Evaluating option IN-OUT:\n")) ;
	        DEBUG2 (("Work->fnzeros "ID" fnpiv "ID" nr_in "ID" nc "ID"\n",
	    	Work->fnzeros, fnpiv, nr_in, nc)) ;
	        DEBUG2 (("fncols "ID" fnrows "ID"\n", fncols, fnrows)) ;

	        /* determine if BLAS-3 update to be applied before extending. */
	        /* update if too many zero entries accumulate in the LU block */
	        fnzeros = Work->fnzeros + fnpiv * (nr_in + nc) ;

	        DEBUG2 (("fnzeros "ID"\n", fnzeros)) ;

	        new_LUsize = (fnpiv+1) * (fnrows + nr_in + fncols + nc) ;

	        DEBUG2 (("new_LUsize "ID"\n", new_LUsize)) ;
	        DEBUG2 (("relax3 %g\n", relax3)) ;

		/* relax3 parameter uses a double relop, but ignore NaN case: */
	        do_update = (((double) fnzeros) / ((double) new_LUsize)) > relax3 ;

	        DEBUG2 (("do_update "ID"\n", do_update))
	    }
	    else
	    {
	    	do_update = TRUE ;
		fnzeros = 0 ;
	        DEBUG2 (("IN-OUT do_update forced true: "ID"\n", do_update))
	    }

	    DEBUG2 (("option IN  OUT: nr "ID" nc "ID" cost "ID"("ID") relax "ID
		"\n", nr_in, nc, thiscost, extra_zeros, do_extend)) ;

	    if (bestcost == EMPTY || thiscost < bestcost)
	    {
		/* this is the best option seen so far */
		Work->pivot_case = IN_OUT ;
		bestcost = thiscost ;
		Work->do_extend = do_extend ;
		Work->do_update = do_update ;
		new_fnzeros = fnzeros ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* construct candidate column not in front, and search for pivot rows */
    /* ---------------------------------------------------------------------- */

    if (bestcost != 0 && pivcol [OUT] != EMPTY)
    {

#ifndef NDEBUG
	DEBUG2 (("out_col column "ID" NOT in front at position = "ID"\n",
		pivcol [OUT], Fcpos [pivcol [OUT]])) ;
	UMF_dump_rowcol (1, Numeric, Work, pivcol [OUT], TRUE) ;
#endif

	/* Find an empty column in the current frontal matrix to use */
	/* as workspace.  There must be one; otherwise, there couldn't be a */
	/* candidate pivot column outside the current front. */

	/* Fd: destination of final pivot column, currently unoccupied */
	/* Use Fd as temporary workspace to construct the pivcol [OUT] */
	Fd = Fx + (fncols_max - fnpiv - 1) * fnrows_max ;

	ASSERT (fncols + fnpiv < fncols_max) ;

	/* ------------------------------------------------------------------ */
	/* construct the candidate column (currently not in the front) */
	/* ------------------------------------------------------------------ */

	/* Construct the column in Fd, Wm, using Wp for the positions: */
	/* Wm [0..cdeg [OUT]-1]	list of row indices in the column */
	/* Fd [0..cdeg [OUT]-1]	list of corresponding numerical values */
	/* Wp [0..n-1] starts as all negative, and ends that way too. */

	cdeg [OUT] = 0 ;

#ifndef NDEBUG
	/* check Wp */
	DEBUG5 (("COL ASSEMBLE: cdeg 0\nREDUCE COL out "ID"\n", pivcol [OUT])) ;
	if (UMF_debug > 0 || MAX (n_row, n_col) < 1000)
	{
	    for (i = 0 ; i < MAX (n_row, n_col) ; i++)
	    {
		ASSERT (Wp [i] < 0) ;
		ASSERT (Wp [i] > Work->Wpflag) ;
	    }
	}
	DEBUG5 (("max_cdeg: "ID"\n", max_cdeg)) ;
#endif

	ASSERT (pivcol [OUT] >= 0 && pivcol [OUT] < n_col) ;
	ASSERT (NON_PIVOTAL_COL (pivcol [OUT])) ;

	tpi = Col_tuples [pivcol [OUT]] ;
	if (tpi)
	{
	    tp = (Tuple *) (Memory + tpi) ;
	    tp1 = tp ;
	    tp2 = tp ;
	    tpend = tp + Col_tlen [pivcol [OUT]] ;
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
		if (Cols [f] == EMPTY) continue ; /* column already assembled */
		ASSERT (pivcol [OUT] == Cols [f]) ;

		Rows = Cols + ep->ncols ;
		nrows = ep->nrows ;
		p += UNITS (Int, ep->ncols + nrows) ;
		C = ((Entry *) p) + f * nrows ;

		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    if (row >= 0) /* skip this if already gone from element */
		    {
			ASSERT (row < n_row) ;
			pos = Wp [row] ;
			if (pos < 0)
			{
			    if (cdeg [OUT] >= max_cdeg)
			    {
				return (UMFPACK_ERROR_different_pattern) ;
			    }
			    /* new entry in the pattern - save Wp */
			    Wp [row] = cdeg [OUT] ;
			    Wm [cdeg [OUT]] = row ;
			    Fd [cdeg [OUT]++] = C [i] ;
			}
			else
			{
			    /* entry already in pattern - sum the values */
			    /* Fd [pos] += C [i] ; */
			    ASSEMBLE (Fd [pos], C [i]) ;
			}
		    }
		}

		*tp2++ = *tp ;	/* leave the tuple in the list */
	    }
	    Col_tlen [pivcol [OUT]] = tp2 - tp1 ;
	}

	/* ------------------------------------------------------------------ */

#ifndef NDEBUG
	DEBUG4 (("Reduced column: cdeg out "ID"\n", cdeg [OUT])) ;
	for (i = 0 ; i < cdeg [OUT] ; i++)
	{
	    DEBUG6 ((" "ID" "ID" "ID, i, Wm [i], Wp [Wm [i]])) ;
	    EDEBUG6 (Fd [i]) ;
	    DEBUG6 (("\n")) ;
	    ASSERT (i == Wp [Wm [i]]) ;
	}
#endif

	/* ------------------------------------------------------------------ */
	/* Clear Wp */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < cdeg [OUT] ; i++)
	{
	    Wp [Wm [i]] = EMPTY ;	/* clear Wp */
	}

	/* ------------------------------------------------------------------ */
	/* new degree of pivcol [OUT] is cdeg [OUT] */
	/* ------------------------------------------------------------------ */

	/* search for two candidate pivot rows */
	status = UMF_row_search (Numeric, Work, Symbolic, cdeg [OUT], Wm,
	    pivrow [OUT], rdeg [OUT], Woi, Woo, pivrow [IN], Fd, pivcol [OUT],
	    freebie) ;

	/* ------------------------------------------------------------------ */
	/* fatal error if matrix pattern has changed since symbolic analysis */
	/* ------------------------------------------------------------------ */

	if (status == UMFPACK_ERROR_different_pattern)
	{
	    return (UMFPACK_ERROR_different_pattern) ;
	}

	/* ------------------------------------------------------------------ */
	/* check for rectangular, singular matrix */
	/* ------------------------------------------------------------------ */

	if (status == UMFPACK_WARNING_singular_matrix)
	{
	    /* Pivot column is empty, and row-merge set is empty too. */
	    /* The matrix is structurally singular.  */

	    DEBUG0 (("Warning: pivcol [OUT]: "ID" discard\n", pivcol [OUT])) ;

	    /* remove the failed pivcol [OUT] from candidate set */
	    jj = jcand [OUT] ;
	    ASSERT (jj >= 0 && jj < jmax) ;
	    ASSERT (pivcol [OUT] == Work->Candidates [jj]) ;
	    if (Work->ncand > MAX_CANDIDATES)
	    {
		Work->Candidates [jj] = Work->nextcand++ ;
	    }
	    else
	    {
		col = Work->Candidates [Work->ncand - 1] ;
		Work->Candidates [jj] = col ;
		Work->Candidates [Work->ncand - 1] = 0 ;
		if (col == pivcol [IN])
		{
		    ASSERT (jcand [IN] == Work->ncand-1) ;
		    jcand [IN] = jj ;
		}
	    }
	    Work->ncand-- ;

	    /* delete all of the tuples, and all contributions to this column */
	    DEBUG1 (("Prune tuples of dead outcol: "ID"\n", pivcol [OUT])) ;
	    Col_tlen [pivcol [OUT]] = 0 ;
	    UMF_mem_free_tail_block (Numeric, Col_tuples [pivcol [OUT]]) ;
	    Col_tuples [pivcol [OUT]] = 0 ;

	    if (pivcol [IN] == EMPTY)
	    {
		/* no pivot found at all */
		return (UMFPACK_WARNING_singular_matrix) ;
	    }
	}

	/* ------------------------------------------------------------------ */

	/* Fd [0 .. cdeg[OUT]-1], the workspace column in the frontal matrix, */
	/* no longer needed */

	if (freebie [IN])
	{
	    /* the "in" row is the same as the "in" row for the "in" column */
	    Woi = Fcols ;
	    rdeg [OUT][IN] = rdeg [IN][IN] ;
	    DEBUG4 (("Freebie in, row "ID"\n", pivrow [IN][IN])) ;
	}

	if (freebie [OUT])
	{
	    /* the "out" row is the same as the "out" row for the "in" column */
	    Woo = Wio ;
	    rdeg [OUT][OUT] = rdeg [IN][OUT] ;
	    DEBUG4 (("Freebie out, row "ID"\n", pivrow [IN][OUT])) ;
	}

	/* ------------------------------------------------------------------ */
	/* evaluate OUT_IN option */
	/* ------------------------------------------------------------------ */

	if (pivrow [OUT][IN] != EMPTY)
	{
	    /* the current front would become an Lson of the new front */

	    did_rowmerge = (cdeg [OUT] == 0) ;
	    if (did_rowmerge)
	    {
		/* pivrow [OUT][IN] was found via row merge search */
		/* it is not (yet) in the pivot column pattern (add it now) */
		if (cdeg [OUT] >= max_cdeg)
		{
		    return (UMFPACK_ERROR_different_pattern) ;
		}
		/* new entry in the pattern */
		Wm [0] = pivrow [OUT][IN] ;
		cdeg [OUT] = 1 ;
		ASSERT (nr_out == EMPTY) ;
	    }

	    nc = rdeg [OUT][IN] - fncols ;
	    ASSERT (nc >= 1) ;

	    /* count rows not in current front */
	    nr_out = 0 ;
#ifndef NDEBUG
	    debug_ok = FALSE ;
#endif
	    for (i = 0 ; i < cdeg [OUT] ; i++)
	    {
		row = Wm [i] ;
		ASSERT (row >= 0 && row < n_row) ;
		ASSERT (NON_PIVOTAL_ROW (row)) ;
		if (Frpos [row] < 0 || Frpos [row] >= fnrows) nr_out++ ;
#ifndef NDEBUG
		/* we must see the pivot row somewhere */
		if (row == pivrow [OUT][IN])
		{
		    ASSERT (Frpos [row] >= 0) ;
		    debug_ok = TRUE ;
		}
#endif
	    }
	    ASSERT (debug_ok) ;

	    thiscost =
		/* each column in front grows by nr_out: */
		(nr_out * fncols) +
		/* new cols not affected by front: */
		((nc-1) * (cdeg [OUT]-1)) ;

	    /* check the cost of relaxed OUT_IN amalgamation */

	    extra_rows = ((fnrows-1) + nr_out) - (cdeg [OUT] - 1) ;
	    ASSERT (extra_rows >= 0) ;
	    ASSERT (fnrows + nr_out == extra_rows + cdeg [OUT]) ;
	    extra_zeros = (nc-1) * extra_rows ;	/* symbolic fill-in */

	    ASSERT (fnrows + nr_out == cdeg [OUT] + extra_rows) ;
	    ASSERT (fncols + nc == rdeg [OUT][IN]) ;

	    /* size of relaxed front: */
	    relaxed_front = (fncols + (nc-1)) * ((fnrows-1) + nr_out) ;

	    /* do relaxed amalgamation if the extra zeros are no more */
	    /* than a fraction (default 0.25) of the relaxed front */
	    /* if relax = 0: no extra zeros allowed */
	    /* if relax = +inf: always amalgamate */
	    if (did_rowmerge)
	    {
		do_extend = FALSE ;
	    }
	    else
	    {
		/* relax parameter uses a double relop, but ignore NaN case: */
		do_extend = ((double) extra_zeros)
		< (relax * (double) relaxed_front) ;
	    }

	    if (do_extend)
	    {
		/* count the cost of relaxed amalgamation */
		thiscost += extra_zeros ;

	        /* determine if BLAS-3 update to be applied before extending. */
	        /* update if too many zero entries accumulate in the LU block */

	        DEBUG2 (("Evaluating option OUT-IN:\n")) ;
	        DEBUG2 ((" Work->fnzeros "ID" fnpiv "ID" nr_out "ID" nc "ID"\n",
	    	Work->fnzeros, fnpiv, nr_out, nc)) ;
	        DEBUG2 (("fncols "ID" fnrows "ID"\n", fncols, fnrows)) ;

	        fnzeros = Work->fnzeros + fnpiv * (nr_out + nc) ;

	        DEBUG2 (("fnzeros "ID"\n", fnzeros)) ;

	        new_LUsize = (fnpiv+1) * (fnrows + nr_out + fncols + nc) ;

	        DEBUG2 (("new_LUsize "ID"\n", new_LUsize)) ;
	        DEBUG2 (("relax3 %g\n", relax3)) ;

		/* relax3 parameter uses a double relop, but ignore NaN case: */
	        do_update = (((double) fnzeros) / ((double) new_LUsize)) > relax3 ;

	        DEBUG2 (("do_update "ID"\n", do_update))

	    }
	    else
	    {
	    	do_update = TRUE ;
		fnzeros = 0 ;
	        DEBUG2 (("OUT-IN do_update forced true: "ID"\n", do_update))
	    }

	    DEBUG2 (("option OUT IN : nr "ID" nc "ID" cost "ID"("ID") relax "ID
		"\n", nr_out, nc, thiscost, extra_zeros, do_extend)) ;

	    if (bestcost == EMPTY || thiscost < bestcost)
	    {
		/* this is the best option seen so far */
		Work->pivot_case = OUT_IN ;
		bestcost = thiscost ;
		Work->do_extend = do_extend ;
		Work->do_update = do_update ;
		new_fnzeros = fnzeros ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* evaluate OUT_OUT option */
	/* ------------------------------------------------------------------ */

	if (pivrow [OUT][OUT] != EMPTY)
	{

	    did_rowmerge = (cdeg [OUT] == 0) ;
	    if (did_rowmerge)
	    {
		/* pivrow [OUT][OUT] was found via row merge search */
		/* it is not (yet) in the pivot column pattern (add it now) */
		if (cdeg [OUT] >= max_cdeg)
		{
		    return (UMFPACK_ERROR_different_pattern) ;
		}
		/* new entry in the pattern */
		Wm [0] = pivrow [OUT][OUT] ;
		cdeg [OUT] = 1 ;
		ASSERT (nr_out == EMPTY) ;
		nr_out = 1 ;
	    }

	    /* count rows not in current front */
	    if (nr_out == EMPTY)
	    {
		nr_out = 0 ;
#ifndef NDEBUG
		debug_ok = FALSE ;
#endif
		for (i = 0 ; i < cdeg [OUT] ; i++)
		{
		    row = Wm [i] ;
		    ASSERT (row >= 0 && row < n_row) ;
		    ASSERT (NON_PIVOTAL_ROW (row)) ;
		    if (Frpos [row] < 0 || Frpos [row] >= fnrows) nr_out++ ;
#ifndef NDEBUG
		    /* we must see the pivot row somewhere */
		    if (row == pivrow [OUT][OUT])
		    {
			ASSERT (Frpos [row] < 0 || Frpos [row] >= fnrows) ;
			debug_ok = TRUE ;
		    }
#endif
		}
		ASSERT (debug_ok) ;
	    }

	    /* count columns not in current front */
	    nc = 0 ;
#ifndef NDEBUG
	    debug_ok = FALSE ;
#endif
	    for (i = 0 ; i < rdeg [OUT][OUT] ; i++)
	    {
		col = Woo [i] ;
		ASSERT (col >= 0 && col < n_col) ;
		ASSERT (NON_PIVOTAL_COL (col)) ;
		if (Fcpos [col] < 0) nc++ ;
#ifndef NDEBUG
		/* we must see the pivot column somewhere */
		if (col == pivcol [OUT])
		{
		    ASSERT (Fcpos [col] < 0) ;
		    debug_ok = TRUE ;
		}
#endif
	    }
	    ASSERT (debug_ok) ;

	    extra_cols = (fncols + (nc-1)) - (rdeg [OUT][OUT] - 1) ;
	    extra_rows = (fnrows + (nr_out-1)) - (cdeg [OUT] - 1) ;
	    ASSERT (extra_rows >= 0) ;
	    ASSERT (extra_cols >= 0) ;
	    extra_zeros = ((nc-1) * extra_rows) + ((nr_out-1) * extra_cols) ;

	    ASSERT (fnrows + nr_out == cdeg [OUT] + extra_rows) ;
	    ASSERT (fncols + nc == rdeg [OUT][OUT] + extra_cols) ;

	    thiscost =
		/* new columns: */
		((nc-1) * (cdeg [OUT]-1)) +
		/* old columns in front grow by nr_out-1: */
		((nr_out-1) * (fncols - extra_cols)) ;

	    /* size of relaxed front: */
	    relaxed_front = (fncols + (nc-1)) * (fnrows + (nr_out-1)) ;

	    /* do relaxed amalgamation if the extra zeros are no more */
	    /* than a fraction (default 0.25) of the relaxed front */
	    /* if relax = 0: no extra zeros allowed */
	    /* if relax = +inf: always amalgamate */
	    if (did_rowmerge)
	    {
	    	do_extend = FALSE ;
	    }
	    else
	    {
		/* relax parameter uses a double relop, but ignore NaN case: */
		do_extend = ((double) extra_zeros)
		< (relax * (double) relaxed_front) ;
	    }

	    if (do_extend)
	    {
		/* count the cost of relaxed amalgamation */
		thiscost += extra_zeros ;

	        DEBUG2 (("Evaluating option OUT-OUT:\n")) ;
	        DEBUG2 (("Work->fnzeros "ID" fnpiv "ID" nr_out "ID" nc "ID"\n",
	    	Work->fnzeros, fnpiv, nr_out, nc)) ;
	        DEBUG2 (("fncols "ID" fnrows "ID"\n", fncols, fnrows)) ;

	        /* determine if BLAS-3 update to be applied before extending. */
	        /* update if too many zero entries accumulate in the LU block */
	        fnzeros = Work->fnzeros + fnpiv * (nr_out + nc) ;

	        DEBUG2 (("fnzeros "ID"\n", fnzeros)) ;

	        new_LUsize = (fnpiv+1) * (fnrows + nr_out + fncols + nc) ;

	        DEBUG2 (("new_LUsize "ID"\n", new_LUsize)) ;
	        DEBUG2 (("relax3 %g\n", relax3)) ;

		/* relax3 parameter uses a double relop, but ignore NaN case: */
	    	do_update = (((double) fnzeros) / ((double) new_LUsize)) > relax3 ;

	        DEBUG2 (("do_update "ID"\n", do_update))

	    }
	    else
	    {
	    	do_update = TRUE ;
		fnzeros = 0 ;
	        DEBUG2 (("OUT-OUT do_update forced true: "ID"\n", do_update))
	    }

	    DEBUG2 (("option OUT OUT: nr "ID" nc "ID" cost "ID"\n",
		rdeg [OUT][OUT], cdeg[OUT], thiscost)) ;

	    if (bestcost == EMPTY || thiscost < bestcost)
	    {
		/* this is the best option seen so far */
		Work->pivot_case = OUT_OUT ;
		bestcost = thiscost ;
		Work->do_extend = do_extend ;
		Work->do_update = do_update ;
		new_fnzeros = fnzeros ;
	    }
	}
    }

    /* At this point, a structural pivot has been found. */
    /* It may be numerically zero, however. */
    ASSERT (Work->pivot_case != EMPTY) ;
    DEBUG2 (("local seach, best option "ID", best cost "ID"\n",
	Work->pivot_case, bestcost)) ;

    /* ====================================================================== */
    /* Pivot row and column, and extension, now determined */
    /* ====================================================================== */

    Work->fnzeros = new_fnzeros ;

    /* ---------------------------------------------------------------------- */
    /* finalize the pivot row and column */
    /* ---------------------------------------------------------------------- */

    switch (Work->pivot_case)
    {
	case IN_IN:
	    DEBUG2 (("IN-IN option selected\n")) ;
	    Work->pivcol_in_front = TRUE ;
	    Work->pivrow_in_front = TRUE ;
	    Work->pivcol = pivcol [IN] ;
	    Work->pivrow = pivrow [IN][IN] ;
	    Work->Wcol = (Int *) NULL ;	/* not accessed */
	    Work->ccdeg = nr_in ;
	    Work->Wrow = Fcols ;
	    Work->rrdeg = rdeg [IN][IN] ;
	    jj = jcand [IN] ;
	    break ;

	case IN_OUT:
	    DEBUG2 (("IN-OUT option selected\n")) ;
	    Work->pivcol_in_front = TRUE ;
	    Work->pivrow_in_front = FALSE ;
	    Work->pivcol = pivcol [IN] ;
	    Work->pivrow = pivrow [IN][OUT] ;
	    Work->Wcol = (Int *) NULL ;	/* not accessed */
	    Work->ccdeg = nr_in ;
	    Work->Wrow = Wio ;
	    Work->rrdeg = rdeg [IN][OUT] ;
	    jj = jcand [IN] ;
	    break ;

	case OUT_IN:
	    DEBUG2 (("OUT-IN option selected\n")) ;
	    Work->pivcol_in_front = FALSE ;
	    Work->pivrow_in_front = TRUE ;
	    Work->pivcol = pivcol [OUT] ;
	    Work->pivrow = pivrow [OUT][IN] ;
	    Work->Wcol = Wm ;
	    Work->ccdeg = cdeg [OUT] ;
	    /* Wrow might be equivalenced to Fcols (Freebie in): */
	    Work->Wrow = Woi ;
	    Work->rrdeg = rdeg [OUT][IN] ;
	    /* Work->Wrow[0..fncols-1] is not there.  See Fcols instead */
	    jj = jcand [OUT] ;
	    break ;

	case OUT_OUT:
	    DEBUG2 (("OUT-OUT option selected\n")) ;
	    Work->pivcol_in_front = FALSE ;
	    Work->pivrow_in_front = FALSE ;
	    Work->pivcol = pivcol [OUT] ;
	    Work->pivrow = pivrow [OUT][OUT] ;
	    Work->Wcol = Wm ;
	    Work->ccdeg = cdeg [OUT] ;
	    /* Wrow might be equivalenced to Wio (Freebie out): */
	    Work->Wrow = Woo ;
	    Work->rrdeg = rdeg [OUT][OUT] ;
	    jj = jcand [OUT] ;
	    break ;

    }

    Wcol = Work->Wcol ;

    if (!Work->pivcol_in_front && pivcol [IN] != EMPTY)
    {
	/* clear Frpos if pivcol [IN] was searched, but not selected */
	for (i = fnrows ; i < cdeg [IN] ; i++)
	{
	    Frpos [Frows [i]] = EMPTY;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* remove pivot column from candidate pivot column list */
    /* ---------------------------------------------------------------------- */

    ASSERT (jj >= 0 && jj < jmax) ;
    ASSERT (Work->pivcol == Work->Candidates [jj]) ;
    if (Work->ncand > MAX_CANDIDATES)
    {
	Work->Candidates [jj] = Work->nextcand++ ;
    }
    else
    {
	Work->Candidates [jj] = Work->Candidates [Work->ncand - 1] ;
	Work->Candidates [Work->ncand - 1] = 0 ;
    }
    Work->ncand-- ;

    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    if (!Work->pivcol_in_front)
    {
	DEBUG2 (("All of Wcol, size "ID"\n", Work->ccdeg)) ;
	ASSERT (Wcol == Work->Wm) ;
	for (i = 0 ; i < Work->ccdeg ; i++)
	{
	    row = Wcol [i] ;
	    DEBUG2 (("Wcol row "ID" Frpos[row] "ID"\n", row, Frpos [row])) ;
	}
    }
    else
    {
	DEBUG2 (("All of old Frows, size "ID"\n", fnrows)) ;
	for (i = 0 ; i < fnrows ; i++)
	{
	    row = Frows [i] ;
	    DEBUG2 (("old Frows row "ID" Frpos[row] "ID"\n", row, Frpos [row]));
	}
	DEBUG2 (("All of (new part of Frows), size "ID"\n",Work->ccdeg));
	for (i = fnrows ; i < fnrows + Work->ccdeg ; i++)
	{
	    row = Frows [i] ;
	    DEBUG2 (("new Frows row "ID" Frpos[row] "ID"\n", row, Frpos [row])) ;
	}
    }
    DEBUG2 (("All of Wrow, size "ID"\n", Work->rrdeg)) ;
    if (Work->pivrow_in_front)
    {
	for (i = 0 ; i < fncols ; i++)
	{
	    col = Fcols [i] ;
	    DEBUG2 (("Wrow col "ID" Fcpos[col] "ID"\n", col, Fcpos [col])) ;
	}
	DEBUG2 (("---\n")) ;
	for (i = fncols ; i < Work->rrdeg ; i++)
	{
	    col = Work->Wrow [i] ;
	    DEBUG2 (("Wrow col "ID" Fcpos[col] "ID"\n", col, Fcpos [col])) ;
	}
    }
    else
    {
	for (i = 0 ; i < Work->rrdeg ; i++)
	{
	    col = Work->Wrow [i] ;
	    DEBUG2 (("Wrow col "ID" Fcpos[col] "ID"\n", col, Fcpos [col])) ;
	}
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* remove pivot row index from pattern of pivot column pattern, */
    /* unless extend_front needs it.  Compress the pivot column pattern. */
    /* ---------------------------------------------------------------------- */

    src = 0 ;
    dest = 0 ;

    if (Work->pivcol_in_front)
    {

	pos = Frpos [Work->pivrow] ;
#ifndef NDEBUG
	if ( Work->pivrow_in_front)
	{
	    ASSERT (pos < fnrows && pos >= 0) ;
	}
	else
	{
	    ASSERT (pos < cdeg [IN] && pos >= fnrows) ;
	}
#endif

	if (!Work->pivrow_in_front)
	{
	    ASSERT (Work->ccdeg > 0) ;
	    ASSERT (cdeg [IN] > fnrows) ;

	    /* remove the pivot row index by shifting the last entry into */
	    /* its position */
	    cdeg [IN]-- ;
	    Work->ccdeg-- ;
	    row = Frows [cdeg [IN]] ;
	    Frows [pos] = row ;
	    Frpos [row] = pos ;

	}

    }
    else
    {
	if (Work->do_extend)
	{
	    /* all cases, when front is being extended */
	    for ( ; src < Work->ccdeg ; src++)
	    {
		ASSERT (Wcol == Work->Wm) ;
		row = Wcol [src] ;
		DEBUG2 (("Pruning Wcol, row "ID" (do_extend)", row)) ;
		if (row != Work->pivrow && Frpos [row] < 0)
		{
		    DEBUG2 ((" keep")) ;
		    Wcol [dest++] = row ;
		}
		DEBUG2 (("\n")) ;
		ASSERT (IMPLIES (Work->pivcol_in_front, Frpos [row] < 0)) ;
	    }
	    /* now Wcol contains only the new rows */
	}
	else
	{
	    for ( ; src < Work->ccdeg ; src++)
	    {
		ASSERT (Wcol == Work->Wm) ;
		row = Wcol [src] ;
		DEBUG2 (("Pruning Wcol, row "ID" (no extend)", row)) ;
		if (row != Work->pivrow)
		{
		    DEBUG2 ((" keep")) ;
		    Wcol [dest++] = row ;
		}
		DEBUG2 (("\n")) ;
	    }
	}
	Work->ccdeg = dest ;
    }


    /* ---------------------------------------------------------------------- */
    /* determine whether to do scan2-row and scan2-col */
    /* ---------------------------------------------------------------------- */

    if (Work->do_extend)
    {
	Work->do_scan2row = (fncols > 0) ;
	Work->do_scan2col = (fnrows > 0) ;
    }
    else
    {
	Work->do_scan2row = (fncols > 0) && Work->pivrow_in_front ;
	Work->do_scan2col = (fnrows > 0) && Work->pivcol_in_front ;
    }

    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    DEBUG2 (("LOCAL SEARCH DONE: pivot column "ID" pivot row: "ID,
	Work->pivcol, Work->pivrow)) ;
    DEBUG2 ((" do_extend: "ID"\n", Work->do_extend)) ;
    DEBUG2 (("do_update: "ID"\n", Work->do_update)) ;
    UMF_dump_rowcol (0, Numeric, Work, Work->pivrow, TRUE) ;
    DEBUG2 (("Pivot Wrow "ID":\n", Work->pivrow)) ;
    if (Work->pivrow_in_front || Work->do_extend)
    {
	for (i = 0 ; i < fncols ; i++)
	{
	    DEBUG3 ((" col:: "ID" \n", Fcols [i])) ;
	}
    }
    for (i = ((Work->pivrow_in_front) ? fncols : 0) ; i < Work->rrdeg ; i++)
    {
	if (col != Work->pivcol && Fcpos [col] < 0)
	{
	    DEBUG3 ((" col:: "ID" (new)\n", Work->Wrow [i])) ;
	}
    }
    UMF_dump_rowcol (1, Numeric, Work, Work->pivcol, TRUE) ;
    DEBUG2 (("Pivot "ID":\n", Work->pivcol)) ;
    if (Work->pivcol_in_front || Work->do_extend)
    {
	for (i = 0 ; i < fnrows ; i++)
	{
	    DEBUG3 ((" row:: "ID" \n", Frows [i])) ;
	}
    }
    if (Work->pivcol_in_front)
    {
	for (i = fnrows ; i < fnrows + Work->ccdeg ; i++)
	{
	    row = Frows [i] ;
	    DEBUG3 ((" row:: "ID" (new, in Frows already)\n", Frows [i])) ;
	    ASSERT (row != Work->pivrow) ;
	}
    }
    else
    {
	for (i = 0 ; i < Work->ccdeg ; i++)
	{
	    ASSERT (Wcol == Work->Wm) ;
	    row = Wcol [i] ;
	    DEBUG3 ((" row:: "ID" (new)\n", row)) ;
	    ASSERT (row != Work->pivrow) ;
	}
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* Pivot row and column have been found */
    /* ---------------------------------------------------------------------- */

    /*
	If pivot column index is in the front:
		pivot column pattern (or just extension) is in
		Work->Wcol [0..Work->ccdeg-1]
	otherwise
		pivot column pattern is in Frows [0.. fnrows + Work->ccdeg-1],
		and Frpos [row] is already set properly.
	In both cases, the pivot row index has been removed from
	the pivot column pattern.

	------------------------------------------------------------------------

	If pivot row index is in the front:
		pivot row pattern (or just extension) is in
		Work->Wrow [fncols..Work->rrdeg-1]
	otherwise
		pivot row pattern is in Work->Wrow [0..Work->rrdeg-1]
	The pivot row pattern has not been pruned.

    */

    return (UMFPACK_OK) ;
}

