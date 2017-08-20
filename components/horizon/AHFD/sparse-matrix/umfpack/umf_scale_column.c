/* ========================================================================== */
/* === UMF_scale_column ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Scale the current pivot column, and log the permutation.
    Store the LU factors.  Called by the kernel.

    Returns TRUE if successful, FALSE if out of memory.
*/

#include "umf_internal.h"
#include "umf_mem_alloc_head_block.h"
#include "umf_mem_free_tail_block.h"
#include "umf_get_memory.h"

/* ========================================================================== */

GLOBAL Int UMF_scale_column
(
    NumericType *Numeric,
    WorkType *Work
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int i, k, k1, fnrows_max, fnrows, fncols, *Frpos, *Fcpos, pos, row, col,
	pivrow, pivcol, *Frows, *Fcols, *Lpattern, *Upattern, *Lpos, *Upos,
	llen, ulen, fncols_max, fnpiv, uilen, lnz, unz, *Row_tuples,
	*Col_tuples, *Rperm, *Cperm, *Lilen, *Uilen, *Lip, *Uip, *Li, *Ui,
	pivcol_position, newLchain, newUchain, pivrow_position, p, size, lip,
	uip, lnzi, lnzx, unzx, lnz2i, lnz2x, unz2i, unz2x, zero_pivot,
	is_nonzero, nan_pivot ;
    Entry *D, x, pivot_value, *Fx, *Fcol, *Frow, *Lval, *Uval ;
    double d ;

#ifndef NDEBUG
    Int *Col_degree, *Row_degree ;
    UMF_allocfail = FALSE ;
    if (UMF_gprob > 0)			    /* a double relop, but ignore NaN case */
    {
	double rrr = ((double) (rand ( ))) / (((double) RAND_MAX) + 1) ;
	DEBUG4 (("Check random %e %e\n", rrr, UMF_gprob)) ;
	UMF_allocfail = rrr < UMF_gprob ;   /* a double relop, but ignore NaN case */
	if (UMF_allocfail)
	{
	    DEBUG1 (("Random garbage collection (scale_column)\n")) ;
	}
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    fnrows = Work->fnrows ;
    fncols = Work->fncols ;

    /* ---------------------------------------------------------------------- */

    Rperm = Numeric->Rperm ;
    Cperm = Numeric->Cperm ;
    Lpos = Numeric->Lpos ;
    Upos = Numeric->Upos ;
    Lilen = Numeric->Lilen ;
    Uilen = Numeric->Uilen ;

    Lip = Numeric->Lip ;
    Uip = Numeric->Uip ;
    D = Numeric->D ;

    /* ---------------------------------------------------------------------- */

    k = Work->npiv++ ;

    Fx = Work->Fx ;
    fnrows_max = Work->fnrows_max ;
    fncols_max = Work->fncols_max ;
    fnpiv = Work->fnpiv ;
    Frpos = Work->Frpos ;
    Fcpos = Work->Fcpos ;
    Frows = Work->Frows ;
    Fcols = Work->Fcols ;
    pivrow = Work->pivrow ;
    pivcol = Work->pivcol ;

    ASSERT (pivrow >= 0 && pivrow < Work->n_row) ;
    ASSERT (pivcol >= 0 && pivcol < Work->n_col) ;
    ASSERT (k < MIN (Work->n_row, Work->n_col)) ;

#ifndef NDEBUG
    Col_degree = Numeric->Cperm ;	/* for NON_PIVOTAL_COL macro */
    Row_degree = Numeric->Rperm ;	/* for NON_PIVOTAL_ROW macro */
    if (k % 1000 == 0) DEBUG0 (("step "ID"\n", k))  ;
#endif

    Row_tuples = Numeric->Uip ;
    Col_tuples = Numeric->Lip ;

    Lpattern = Work->Lpattern ;
    llen = Work->llen ;
    Upattern = Work->Upattern ;
    ulen = Work->ulen ;

    /* ---------------------------------------------------------------------- */

    /* Frpos [row] >= 0 for each row in pivot column pattern.   */
    /* offset into pattern is given by:			   	*/
    /* Frpos [row] == offset - 1				*/
    /* Frpos [pivrow] is the offset of the latest pivot row	*/

    /* Fcpos [col] >= 0 for each col in pivot row pattern.	*/
    /* Fcpos [col] == (offset - 1) * fnrows_max		 	*/
    /* Fcpos [pivcol] is the offset of the latest pivot column  */

    /* Fcols [0..fncols-1] is the pivot row pattern (excl pivot cols) */
    /* Frows [0..fnrows-1] is the pivot col pattern (excl pivot rows) */

#ifndef NDEBUG
    DEBUG7 (("Current frontal matrix: (prior to pivcol scale)\n")) ;
    UMF_dump_current_front (Numeric, Work, TRUE) ;
#endif

    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    {
	Int row, col, i, lcnt, ucnt ;

	DEBUG2 (("Store column of L, k = "ID", llen "ID"\n", k, llen)) ;
	for (i = 0 ; i < llen ; i++)
	{
	    row = Lpattern [i] ;
	    ASSERT (row >= 0 && row < Work->n_row) ;
	    DEBUG2 (("    Lpattern["ID"] "ID" Lpos "ID, i, row, Lpos [row])) ;
	    if (row == pivrow) DEBUG2 ((" <- pivot row")) ;
	    DEBUG2 (("\n")) ;
	    ASSERT (NON_PIVOTAL_ROW (row)) ;
	    ASSERT (i == Lpos [row]) ;
	}

	DEBUG2 (("Store row of U, k = "ID", ulen "ID"\n", k, ulen)) ;
	for (i = 0 ; i < ulen ; i++)
	{
	    col = Upattern [i] ;
	    DEBUG2 (("    Upattern["ID"] "ID, i, col)) ;
	    if (col == pivcol) DEBUG2 ((" <- pivot col")) ;
	    DEBUG2 (("\n")) ;
	    ASSERT (col >= 0 && col < Work->n_col) ;
	    ASSERT (NON_PIVOTAL_COL (col)) ;
	    ASSERT (i == Upos [col]) ;
	}

	lcnt = 0 ;
	ucnt = 0 ;
	if (Work->n_row < 1000)
	{
	    for (row = 0 ; row < Work->n_row ; row++)
	    {
		if (NON_PIVOTAL_ROW (row) && Lpos [row] != EMPTY) lcnt++ ;
	    }
	    ASSERT (lcnt == llen) ;
	}
	if (Work->n_col < 1000)
	{
	    for (col = 0 ; col < Work->n_col ; col++)
	    {
		if (NON_PIVOTAL_COL (col) && Upos [col] != EMPTY) ucnt++ ;
	    }
	    ASSERT (ucnt == ulen) ;
	}
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* remove pivot row from L */
    /* ---------------------------------------------------------------------- */

    /* remove pivot row index from current column of L */
    /* if a new Lchain starts, then all entries are removed later */
    DEBUG2 (("Removing pivrow from Lpattern, k = "ID"\n", k)) ;
    ASSERT (NON_PIVOTAL_ROW (pivrow)) ;
    pivrow_position = Lpos [pivrow] ;
    if (pivrow_position != EMPTY)
    {
	/* place the last entry in the column in the */
	/* position of the pivot row index */
	ASSERT (pivrow == Lpattern [pivrow_position]) ;
	row = Lpattern [--llen] ;
	ASSERT (NON_PIVOTAL_ROW (row)) ;
	Lpattern [pivrow_position] = row ;
	Lpos [row] = pivrow_position ;
	Lpos [pivrow] = EMPTY ;
    }

    /* ---------------------------------------------------------------------- */
    /* store the pivot value, for the diagonal matrix D */
    /* ---------------------------------------------------------------------- */

    /* fnpiv-th pivot in frontal matrix located in */
    /* Fx (fnrows_max-fnpiv, fncols_max-fnpiv) */

    Fcol = Fx + (fncols_max - fnpiv) * fnrows_max ;
    pivot_value = Fcol [fnrows_max - fnpiv] ;
    D [k] = pivot_value ;

    ABS (d, pivot_value) ;
    zero_pivot = SCALAR_IS_ZERO (d) ;
    nan_pivot = SCALAR_IS_NAN (d) ;

    if (k == 0)
    {
	Numeric->min_udiag = d ;
	Numeric->max_udiag = d ;
    }
    else
    {
	/* min (abs (diag (U))) behaves as follows:  If any entry is zero,
	   then the result is zero (regardless of the presence of NaN's).
	   Otherwise, if any entry is NaN, then the result is NaN.  Otherwise,
	   the result is the smallest absolute value on the diagonal of U.
	*/

	if (SCALAR_IS_NONZERO (Numeric->min_udiag))
	{
	    if (zero_pivot || nan_pivot)
	    {
		Numeric->min_udiag = d ;
	    }
	    else if (!SCALAR_IS_NAN (Numeric->min_udiag))
	    {
		/* d and min_udiag are both non-NaN */
		Numeric->min_udiag = MIN (Numeric->min_udiag, d) ;
	    }
	}

	/*
	   max (abs (diag (U))) behaves as follows:  If any entry is NaN
	   then the result is NaN.  Otherise, the result is the largest
	   absolute value on the diagonal of U.
	*/

	if (nan_pivot)
	{
	    Numeric->max_udiag = d ;
	}
	else if (!SCALAR_IS_NAN (Numeric->max_udiag))
	{
	    /* d and max_udiag are both non-NaN */
	    Numeric->max_udiag = MAX (Numeric->max_udiag, d) ;
	}
    }

    if (!zero_pivot)
    {
	/* the pivot is nonzero, but might be Inf or NaN */
	Numeric->nnzpiv++ ;
    }
    DEBUG4 (("Pivot abs value: %g nnzpiv: "ID" D["ID"]=", d, Numeric->nnzpiv, k)) ;
    EDEBUG4 (pivot_value) ;
    DEBUG4 (("\n")) ;

    /* ---------------------------------------------------------------------- */
    /* scale pivot column and count nonzeros in kth column of L */
    /* ---------------------------------------------------------------------- */

    lnz = 0 ;
    lnz2i = 0 ;
    lnz2x = llen ;
    for (i = 0 ; i < fnrows ; i++)
    {
	EDEBUG4 (Fcol [i]) ;
	EDEBUG4 (pivot_value) ;
	if (IS_NONZERO (Fcol [i]))
	{
	    /* Fcol [i] is nonzero, NaN, or Inf */
	    /* Fcol [i] = Fcol [i] / pivot_value ; */
	    DIV (x, Fcol [i], pivot_value) ;
	    Fcol [i] = x ;
	    /* underflow may have occured, so check if result is zero */
	    is_nonzero = IS_NONZERO (x) ;
	}
	else
	{
	    /* Fcol [i] is zero.  Do not divide by pivot value. */
	    is_nonzero = FALSE ;
	}
	DEBUG4 (("pivot column: "ID" is_nonzero "ID, i, is_nonzero)) ;
	EDEBUG4 (Fcol [i]) ;
	DEBUG4 (("\n")) ;

	/* if we start a new Lchain: */
	if (is_nonzero)
	{
	    /* one new integer and one new Entry */
	    DEBUG4 ((" got an entry \n")) ;
	    lnz++ ;
	}

	/* if we continue the prior Lchain: */
	row = Frows [i] ;
	ASSERT (row != pivrow) ;
	ASSERT (NON_PIVOTAL_ROW (row)) ;
	pos = Lpos [row] ;
	if (pos == EMPTY && is_nonzero)
	{
	    /* row is not in the Lpattern, add it if Fcol [i] is nonzero */
	    lnz2i++ ;
	    lnz2x++ ;
	}
	DEBUG3 (("Scale L col, row "ID" pos "ID" scaled value", row, pos)) ;
	EDEBUG3 (Fcol [i]) ;
	DEBUG3 ((" ::: lnz "ID" lnz2i "ID" lnz2x "ID"\n", lnz, lnz2i, lnz2x)) ;
    }

    /* determine if we start a new Lchain or continue the old one */
    if (llen == 0 || zero_pivot)
    {
	/* llen == 0 means there is no prior Lchain */
	/* D [k] == 0 means the pivot column is empty */
	newLchain = TRUE ;
    }
    else
    {
	newLchain =
		/* storage for starting a new Lchain */
		UNITS (Entry, lnz) + UNITS (Int, lnz)
	    <=
		/* storage for continuing a prior Lchain */
		UNITS (Entry, lnz2x) + UNITS (Int, lnz2i) ;
    }

    if (newLchain)
    {
	/* start a new chain for column k of L */
	DEBUG2 (("Start new Lchain, k = "ID"\n", k)) ;

	pivrow_position = EMPTY ;

	/* clear the prior Lpattern */
	for (i = 0 ; i < llen ; i++)
	{
	    row = Lpattern [i] ;
	    ASSERT (NON_PIVOTAL_ROW (row)) ;
	    Lpos [row] = EMPTY ;
	}
	llen = 0 ;

	lnzi = lnz ;
	lnzx = lnz ;
    }
    else
    {
	/* continue the prior Lchain */
	DEBUG2 (("Continue  Lchain, k = "ID"\n", k)) ;
	lnzi = lnz2i ;
	lnzx = lnz2x ;
    }

    /* ---------------------------------------------------------------------- */
    /* count the nonzeros in the row of U */
    /* ---------------------------------------------------------------------- */

    /* store the numerical entries and find new nonzeros */
    Frow = Fx + (fnrows_max - fnpiv) ;

    unz = 0 ;
    unz2i = 0 ;
    unz2x = ulen ;
    DEBUG2 (("unz2x is "ID"\n", unz2x)) ;

    /* if row k does not end a Uchain, pivcol will not be included in ulen */

    ASSERT (NON_PIVOTAL_COL (pivcol)) ;
    pivcol_position = Upos [pivcol] ;
    if (pivcol_position != EMPTY)
    {
	unz2x-- ;
	DEBUG2 (("(exclude pivcol) unz2x is now "ID"\n", unz2x)) ;
    }

    ASSERT (unz2x >= 0) ;

    for (i = 0 ; i < fncols ; i++)
    {
	x = Frow [i * fnrows_max] ;
	is_nonzero = IS_NONZERO (x) ;

	/* if we start a new Uchain */
	if (is_nonzero)
	{
	    unz++ ;
	    DEBUG2 (("If Unew: "ID, unz)) ;
	    EDEBUG2 (x) ;
	    DEBUG2 (("\n")) ;
	}

	/* if we continue the prior Uchain */
	col = Fcols [i] ;
	ASSERT (col != pivcol) ;
	ASSERT (NON_PIVOTAL_COL (col)) ;
	pos = Upos [col] ;
	if (pos == EMPTY && is_nonzero)
	{
	    /* add this new nonzero entry to the U pattern, if nonzero */
	    unz2i++ ;
	    unz2x++ ;
	    DEBUG2 (("If old:                      "ID" :", unz2x)) ;
	    EDEBUG2 (x) ;
	    DEBUG2 (("\n")) ;
	}
    }

    ASSERT (IMPLIES (k == 0, ulen == 0)) ;

    /* determine if we start a new Uchain or continue the old one */
    if (ulen == 0 || zero_pivot)
    {
	/* ulen == 0 means there is no prior Uchain */
	/* D [k] == 0 means the matrix is singular (pivot row might */
	/* not be empty, however, but start a new Uchain to prune zero */
	/* entries for the deg > 0 test in UMF_u*solve) */
	newUchain = TRUE ;
    }
    else
    {
	newUchain =
		/* approximate storage for starting a new Uchain */
		UNITS (Entry, unz) + UNITS (Int, unz)
	    <=
		/* approximate storage for continuing a prior Uchain */
		UNITS (Entry, unz2x) + UNITS (Int, unz2i) ;

	/* this would be exact, except for the Int to Unit rounding, */
	/* because the Upattern is stored only at the end of the Uchain */
    }

    /* ---------------------------------------------------------------------- */
    /* allocate space for the column of L and the row of U */
    /* ---------------------------------------------------------------------- */

    size = UNITS (Int, lnzi) + UNITS (Entry, lnzx) ;
    if (newUchain)
    {
	/* store the pattern of the last row in the prior Uchain */
	size += UNITS (Int, ulen) ;
	unzx = unz ;
    }
    else
    {
	unzx = unz2x ;
    }
    size += UNITS (Entry, unzx) ;

    p = UMF_mem_alloc_head_block (Numeric, size) ;
    if (!p)
    {
	/* do garbage collection, realloc, and try again */
	if (!UMF_get_memory (Numeric, Work, size))
	{
	    return (FALSE) ;	/* out of memory */
	}
	p = UMF_mem_alloc_head_block (Numeric, size) ;
    }
    if (!p)
    {
	return (FALSE) ;	/* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* store the column of L */
    /* ---------------------------------------------------------------------- */

    lip = p ;

    Li = (Int *) (Numeric->Memory + p) ;
    p += UNITS (Int, lnzi) ;
    Lval = (Entry *) (Numeric->Memory + p) ;
    p += UNITS (Entry, lnzx) ;

    for (i = 0 ; i < lnzx ; i++)
    {
	CLEAR (Lval [i]) ;
    }

    /* store the numerical entries */

    if (newLchain)
    {
	/* flag the first column in the Lchain by negating Lip [k] */
	lip = -lip ;

	ASSERT (llen == 0) ;
	for (i = 0 ; i < fnrows ; i++)
	{
	    x = Fcol [i] ;
	    DEBUG4 (("Store column, i "ID" is_nonzero(x) "ID"\n", i,
		IS_NONZERO (x))) ;
	    EDEBUG4 (x) ;
	    EDEBUG4 (Fcol [i]) ;
	    DEBUG4 (("\n")) ;

	    if (IS_NONZERO (x))
	    {
		DEBUG4 (("Store column, i "ID" is_nonzero(x) "ID"\n", i,
		IS_NONZERO (x))) ;

		row = Frows [i] ;
		ASSERT (NON_PIVOTAL_ROW (row)) ;
		pos = llen++ ;
		Lpattern [pos] = row ;
		Lpos [row] = pos ;
		Li [pos] = row ;
		Lval [pos] = x ;
		DEBUG2 (("(newLchain) New entry in Lpattern: row "ID" pos "ID
		    "\n", row, pos)) ;
		DEBUG2 (("(newLchain) Store Lval row "ID" pos "ID" value",
		    row, pos)) ;
		EDEBUG2 (x) ;
		DEBUG2 (("\n")) ;
	    }
	}
    }
    else
    {
	ASSERT (llen > 0) ;
	for (i = 0 ; i < fnrows ; i++)
	{
	    x = Fcol [i] ;
	    if (IS_NONZERO (x))
	    {
		row = Frows [i] ;
		ASSERT (NON_PIVOTAL_ROW (row)) ;
		pos = Lpos [row] ;
		if (pos == EMPTY)
		{
		    /* add this new nonzero entry to the L pattern */
		    pos = llen++ ;
		    DEBUG2 (("New entry in Lpattern: row "ID" pos "ID"\n",
			row, pos)) ;
		    ASSERT (llen <= lnzx) ;
		    Lpattern [pos] = row ;
		    Lpos [row] = pos ;
		    *Li++ = row ;
		}
		DEBUG2 (("Store Lval row "ID" pos "ID" value", row, pos)) ;
		EDEBUG2 (x) ;
		DEBUG2 (("\n")) ;
		ASSERT (row == Lpattern [pos]) ;
		ASSERT (pos < lnzx) ;
		Lval [pos] = x ;
	    }
	}
    }
    DEBUG4 (("llen "ID" lnzx "ID"\n", llen, lnzx)) ;
    ASSERT (llen == lnzx) ;
    ASSERT (lnz <= llen) ;

    Numeric->lnz += lnz ;
    Numeric->nLentries += lnzx ;
    Work->llen = llen ;
    Numeric->isize += lnzi ;

    /* ---------------------------------------------------------------------- */
    /* store the row of U */
    /* ---------------------------------------------------------------------- */

    uip = p ;

    if (newUchain)
    {
	/* starting a new Uchain - flag this by negating Uip [k] */
	uip = -uip ;
	DEBUG2 (("Start new Uchain, k = "ID"\n", k)) ;

	pivcol_position = EMPTY ;

	/* end the prior Uchain */
	/* save the current Upattern, and then */
	/* clear it and start a new Upattern */
	DEBUG2 (("Ending prior chain, k-1 = "ID"\n", k-1)) ;
	uilen = ulen ;
	Ui = (Int *) (Numeric->Memory + p) ;
	Numeric->isize += ulen ;
	p += UNITS (Int, ulen) ;
	for (i = 0 ; i < ulen ; i++)
	{
	    col = Upattern [i] ;
	    ASSERT (col >= 0 && col < Work->n_col) ;
	    ASSERT (NON_PIVOTAL_COL (col)) ;
	    Upos [col] = EMPTY ;
	    Ui [i] = col ;
	}

	ulen = 0 ;

    }
    else
    {
	/* continue the prior Uchain */
	DEBUG2 (("Continue  Uchain, k = "ID"\n", k)) ;
	ASSERT (k > 0) ;

	/* remove pivot col index from current row of U */
	/* if a new Uchain starts, then all entries are removed later */
	DEBUG2 (("Removing pivcol from Upattern, k = "ID"\n", k)) ;

	if (pivcol_position != EMPTY)
	{
	    /* place the last entry in the row in the */
	    /* position of the pivot col index */
	    ASSERT (pivcol == Upattern [pivcol_position]) ;
	    col = Upattern [--ulen] ;
	    ASSERT (col >= 0 && col < Work->n_col) ;
	    ASSERT (NON_PIVOTAL_COL (col)) ;
	    Upattern [pivcol_position] = col ;
	    Upos [col] = pivcol_position ;
	    Upos [pivcol] = EMPTY ;
	}

	/* this row continues the Uchain.  Keep track of how much */
	/* to trim from the k-th length to get the length of the */
	/* (k-1)st row of U */
	uilen = unz2i ;

    }

    Uval = (Entry *) (Numeric->Memory + p) ;
    /* p += UNITS (Entry, unzx), no need to increment p */

    for (i = 0 ; i < unzx ; i++)
    {
	CLEAR (Uval [i]) ;
    }

    if (newUchain)
    {
	ASSERT (ulen == 0) ;
	for (i = 0 ; i < fncols ; i++)
	{
	    x = Frow [i * fnrows_max] ;
	    if (IS_NONZERO (x))
	    {
		/* add this new nonzero entry to the U pattern */
		col = Fcols [i] ;
		ASSERT (col >= 0 && col < Work->n_col) ;
		ASSERT (NON_PIVOTAL_COL (col)) ;
		pos = ulen++ ;
		Upattern [pos] = col ;
		Upos [col] = pos ;
		Uval [pos] = x ;
		DEBUG2 (("(newUchain) New entry in Upattern: col "ID" pos "ID
		    "\n", col, pos)) ;
		DEBUG2 (("(newUchain) Store Uval col "ID" pos "ID" value",
		    col, pos)) ;
		EDEBUG2 (x) ;
		DEBUG2 (("\n")) ;
	    }
	}
    }
    else
    {

	ASSERT (ulen > 0) ;

	/* store the numerical entries and find new nonzeros */

	for (i = 0 ; i < fncols ; i++)
	{
	    x = Frow [i * fnrows_max] ;
	    if (IS_NONZERO (x))
	    {
		col = Fcols [i] ;
		ASSERT (col >= 0 && col < Work->n_col) ;
		ASSERT (NON_PIVOTAL_COL (col)) ;
		pos = Upos [col] ;
		if (pos == EMPTY)
		{
		    /* add this new nonzero entry to the U pattern */
		    ASSERT (ulen < unzx) ;
		    pos = ulen++ ;
		    Upattern [pos] = col ;
		    Upos [col] = pos ;
		    DEBUG2 (("New entry in Upattern: col "ID" pos "ID"\n",
			col, pos)) ;
		}
		DEBUG2 (("Store Uval col "ID" pos "ID" value", col, pos)) ;
		EDEBUG2 (x) ;
		DEBUG2 (("\n")) ;
		ASSERT (col == Upattern [pos]) ;
		ASSERT (pos < unzx) ;
		Uval [pos] = x ;
	    }
	}
    }

    ASSERT (ulen == unzx) ;
    ASSERT (unz <= ulen) ;
    Numeric->unz += unz ;
    Numeric->nUentries += unzx ;
    Work->ulen = ulen ;
    DEBUG1 (("Work->ulen = "ID" at end of pivot step, k: "ID"\n", ulen, k)) ;

    /* ---------------------------------------------------------------------- */
    /* count the "true" flops, based on LU pattern only */
    /* ---------------------------------------------------------------------- */

    /* the outer product is not done here, but in the BLAS */
    Numeric->flops += DIV_FLOPS * lnz		/* scale pivot column */
	+ MULTSUB_FLOPS * (lnz*unz) ;		/* outer product */

    /* ====================================================================== */
    /* A pivot step is complete */
    /* ====================================================================== */

#ifndef NDEBUG
    DEBUG7 (("Current frontal matrix: (after pivcol scale)\n")) ;
    UMF_dump_current_front (Numeric, Work, TRUE) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* remove pivot row and column from frontal pattern */
    /* ---------------------------------------------------------------------- */

    Frpos [pivrow] = EMPTY ;
    Fcpos [pivcol] = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* deallocate the pivot row and pivot column tuples */
    /* ---------------------------------------------------------------------- */

    UMF_mem_free_tail_block (Numeric, Row_tuples [pivrow]) ;
    UMF_mem_free_tail_block (Numeric, Col_tuples [pivcol]) ;

    /* ---------------------------------------------------------------------- */
    /* the pivot column is fully assembled and scaled, and is now the */
    /* k-th column of L. The pivot row is the k-th row of U. */
    /* ---------------------------------------------------------------------- */

    DEBUG5 (("number of pivots prior to this one: "ID"\n", k)) ;
    ASSERT (NON_PIVOTAL_ROW (pivrow)) ;
    ASSERT (NON_PIVOTAL_COL (pivcol)) ;

    /* save row and column inverse permutation */
    k1 = ONES_COMPLEMENT (k) ;
    Rperm [pivrow] = k1 ;			/* aliased with Row_degree */
    Cperm [pivcol] = k1 ;			/* aliased with Col_degree */

    ASSERT (!NON_PIVOTAL_ROW (pivrow)) ;
    ASSERT (!NON_PIVOTAL_COL (pivcol)) ;

    Lpos [pivrow] = pivrow_position ;
    Upos [pivcol] = pivcol_position ;

    Lip [pivcol] = lip ;			/* aliased with Col_tuples */
    Lilen [pivcol] = lnzi ;			/* aliased with Col_tlen */

    Uip [pivrow] = uip ;			/* aliased with Row_tuples */
    Uilen [pivrow] = uilen ;			/* aliased with Row_tlen */

    return (TRUE) ;

}

