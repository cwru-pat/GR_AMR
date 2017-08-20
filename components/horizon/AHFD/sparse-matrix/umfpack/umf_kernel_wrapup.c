/* ========================================================================== */
/* === UMF_kernel_wrapup ==================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* The matrix is factorized.  Finish the LU data structure. */

#include "umf_internal.h"

GLOBAL void UMF_kernel_wrapup
(
    NumericType *Numeric,
    SymbolicType *Symbolic,
    WorkType *Work
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int i, k, col, row, llen, ulen, *ip, *Rperm, *Cperm, *Lilen, npiv, lp,
	*Uilen, *Lip, *Uip, *Cperm_init, up, pivrow, pivcol, *Lpos, *Upos, *Wr,
	*Wc, *Wp, *Frpos, *Fcpos, *Row_degree, *Col_degree, *Rperm_init,
	n_row, n_col, n_inner ;
    Entry *D ;

#ifndef NDEBUG
    UMF_dump_matrix (Numeric, Work, FALSE) ;
#endif

    DEBUG0 (("Kernel complete, Starting Kernel wrapup\n")) ;
    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;
    n_inner = MIN (n_row, n_col) ;
    Rperm = Numeric->Rperm ;
    Cperm = Numeric->Cperm ;
    Lilen = Numeric->Lilen ;
    Uilen = Numeric->Uilen ;
    Upos = Numeric->Upos ;
    Lpos = Numeric->Lpos ;
    Lip = Numeric->Lip ;
    Uip = Numeric->Uip ;
    D = Numeric->D ;

    npiv = Work->npiv ;
    Numeric->npiv = npiv ;
    Numeric->ulen = Work->ulen ;

    ASSERT (n_row == Numeric->n_row) ;
    ASSERT (n_col == Symbolic->n_col) ;
    DEBUG0 (("Wrap-up: npiv "ID" ulen "ID"\n", npiv, Numeric->ulen)) ;
    ASSERT (npiv <= n_inner) ;

    /* this will be nonzero only if matrix is singular or rectangular */
    ASSERT (IMPLIES (npiv == n_col, Work->ulen == 0)) ;

    /* ---------------------------------------------------------------------- */
    /* check if matrix is singular or rectangular */
    /* ---------------------------------------------------------------------- */

    Col_degree = Cperm ;	/* for NON_PIVOTAL_COL macro */
    Row_degree = Rperm ;	/* for NON_PIVOTAL_ROW macro */

    if (npiv < n_row)
    {
	/* finalize the row permutation */
	k = npiv ;
	for (row = 0 ; row < n_row ; row++)
	{
	    if (NON_PIVOTAL_ROW (row))
	    {
		Rperm [row] = ONES_COMPLEMENT (k) ;
		ASSERT (!NON_PIVOTAL_ROW (row)) ;
		Lpos [row] = EMPTY ;
		Uip [row] = EMPTY ;
		Uilen [row] = 0 ;
		k++ ;
	    }
	}
	ASSERT (k == n_row) ;
    }

    if (npiv < n_col)
    {
	/* finalize the col permutation */
	k = npiv ;
	for (col = 0 ; col < n_col ; col++)
	{
	    if (NON_PIVOTAL_COL (col))
	    {
		Cperm [col] = ONES_COMPLEMENT (k) ;
		ASSERT (!NON_PIVOTAL_COL (col)) ;
		Upos [col] = EMPTY ;
		Lip [col] = EMPTY ;
		Lilen [col] = 0 ;
		k++ ;
	    }
	}
	ASSERT (k == n_col) ;
    }

    if (npiv < n_inner)
    {
	/* finalize the diagonal of U */
	for (k = npiv ; k < n_inner ; k++)
	{
	    CLEAR (D [k]) ;
	}
    }

    /* save the pattern of the last row of U */
    if (Numeric->ulen > 0)
    {
	Numeric->Upattern = Work->Upattern ;
	Work->Upattern = (Int *) NULL ;
    }

    ASSERT (Numeric->nnzpiv <= npiv) ;
    if (Numeric->nnzpiv < n_inner)
    {
	/* the rest of the diagonal is zero, so min_udiag becomes 0 */
	Numeric->min_udiag = 0.0 ;
    }

    /* ---------------------------------------------------------------------- */
    /* size n_row, n_col workspaces that can be used here: */
    /* ---------------------------------------------------------------------- */

    Frpos = Work->Frpos ;	/* of size n_row+1 */
    Fcpos = Work->Fcpos ;	/* of size n_col+1 */
    Wp = Work->Wp ;		/* of size MAX(n_row,n_col)+1 */
    /* Work->Upattern ;		cannot be used (in Numeric) */
    Wr = Work->Lpattern ;	/* of size n_row+1 */
    Wc = Work->E ;		/* of size n_col+n_inner+1 */

    /* ---------------------------------------------------------------------- */
    /* construct Rperm from inverse permutations */
    /* ---------------------------------------------------------------------- */

    /* use Frpos for temporary copy of inverse row permutation [ */

    for (pivrow = 0 ; pivrow < n_row ; pivrow++)
    {
	k = Rperm [pivrow] ;
	ASSERT (k < 0) ;
	k = ONES_COMPLEMENT (k) ;
	ASSERT (k >= 0 && k < n_row) ;
	Wp [k] = pivrow ;
	Frpos [pivrow] = k ;
    }
    for (k = 0 ; k < n_row ; k++)
    {
	Rperm [k] = Wp [k] ;
    }

    /* ---------------------------------------------------------------------- */
    /* construct Cperm from inverse permutation */
    /* ---------------------------------------------------------------------- */

    /* use Fcpos for temporary copy of inverse column permutation [ */

    for (pivcol = 0 ; pivcol < n_col ; pivcol++)
    {
	k = Cperm [pivcol] ;
	ASSERT (k < 0) ;
	k = ONES_COMPLEMENT (k) ;
	ASSERT (k >= 0 && k < n_col) ;
	Wp [k] = pivcol ;
	/* save a copy of the inverse column permutation in Fcpos */
	Fcpos [pivcol] = k ;
    }
    for (k = 0 ; k < n_col ; k++)
    {
	Cperm [k] = Wp [k] ;
    }

#ifndef NDEBUG
    for (k = 0 ; k < n_col ; k++)
    {
	col = Cperm [k] ;
	ASSERT (col >= 0 && col < n_col) ;
	ASSERT (Fcpos [col] == k) ;		/* col is the kth pivot */
    }
    for (k = 0 ; k < n_row ; k++)
    {
	row = Rperm [k] ;
	ASSERT (row >= 0 && row < n_row) ;
	ASSERT (Frpos [row] == k) ;		/* row is the kth pivot */
    }
#endif

#ifndef NDEBUG
    UMF_dump_lu (Numeric) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* permute Lpos, Upos, Lilen, Lip, Uilen, and Uip */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < npiv ; k++)
    {
	pivrow = Rperm [k] ;
	Wr [k] = Uilen [pivrow] ;
	Wp [k] = Uip [pivrow] ;
    }

    for (k = 0 ; k < npiv ; k++)
    {
	Uilen [k] = Wr [k] ;
	Uip [k] = Wp [k] ;
    }

    for (k = 0 ; k < npiv ; k++)
    {
	pivrow = Rperm [k] ;
	Wp [k] = Lpos [pivrow] ;
    }

    for (k = 0 ; k < npiv ; k++)
    {
	Lpos [k] = Wp [k] ;
    }

    for (k = 0 ; k < npiv ; k++)
    {
	pivcol = Cperm [k] ;
	Wc [k] = Lilen [pivcol] ;
	Wp [k] = Lip [pivcol] ;
    }

    for (k = 0 ; k < npiv ; k++)
    {
	Lilen [k] = Wc [k] ;
	Lip [k] = Wp [k] ;
    }

    for (k = 0 ; k < npiv ; k++)
    {
	pivcol = Cperm [k] ;
	Wp [k] = Upos [pivcol] ;
    }

    for (k = 0 ; k < npiv ; k++)
    {
	Upos [k] = Wp [k] ;
    }

    /* ---------------------------------------------------------------------- */
    /* terminate the last Uchain and last Lchain */
    /* ---------------------------------------------------------------------- */

    Upos [npiv] = EMPTY ;
    Lpos [npiv] = EMPTY ;
    Uip [npiv] = EMPTY ;
    Lip [npiv] = EMPTY ;
    Uilen [npiv] = 0 ;
    Lilen [npiv] = 0 ;

    /* ---------------------------------------------------------------------- */
    /* convert U to the new pivot order */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < npiv ; k++)
    {
	up = Uip [k] ;
	if (up < 0)
	{
	    /* this is the start of a new Uchain (with a pattern) */
	    ulen = Uilen [k] ;
	    DEBUG4 (("K "ID" New U.  ulen "ID" End_Uchain 1\n", k, ulen)) ;
	    if (ulen > 0)
	    {
		up = -up ;
		ip = (Int *) (Numeric->Memory + up) ;
		for (i = 0 ; i < ulen ; i++)
		{
		    col = *ip ;
		    DEBUG4 (("    old col "ID" new col "ID"\n", col, Fcpos [col]));
		    ASSERT (col >= 0 && col < n_col) ;
		    *ip++ = Fcpos [col] ;
		}
	    }
	}
    }

    ulen = Numeric->ulen ;
    if (ulen > 0)
    {
	/* convert last pivot row of U to the new pivot order */
	DEBUG4 (("K "ID" (last)\n")) ;
	for (i = 0 ; i < ulen ; i++)
	{
	    col = Numeric->Upattern [i] ;
	    DEBUG4 (("    old col "ID" new col "ID"\n", col, Fcpos [col])) ;
	    Numeric->Upattern [i] = Fcpos [col] ;
	}
    }

    /* Fcpos no longer needed ] */

    /* ---------------------------------------------------------------------- */
    /* convert L to the new pivot order */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < npiv ; k++)
    {
	llen = Lilen [k] ;
	DEBUG4 (("K "ID" New L.  llen "ID" \n", k, llen)) ;
	if (llen > 0)
	{
	    lp = Lip [k] ;
	    if (lp < 0)
	    {
	    	/* this starts a new Lchain */
		lp = -lp ;
	    }
	    ip = (Int *) (Numeric->Memory + lp) ;
	    for (i = 0 ; i < llen ; i++)
	    {
		row = *ip ;
		DEBUG4 (("    old row "ID" new row "ID"\n", row, Frpos [row])) ;
		ASSERT (row >= 0 && row < n_row) ;
		*ip++ = Frpos [row] ;
	    }
	}
    }

    /* Frpos no longer needed ] */

    /* ---------------------------------------------------------------------- */
    /* combine symbolic and numeric permutations */
    /* ---------------------------------------------------------------------- */

    Cperm_init = Symbolic->Cperm_init ;
    Rperm_init = Symbolic->Rperm_init ;

    for (k = 0 ; k < n_row ; k++)
    {
	Rperm [k] = Rperm_init [Rperm [k]] ;
    }

    for (k = 0 ; k < n_col ; k++)
    {
        Cperm [k] = Cperm_init [Cperm [k]] ;
    }

    /* Work object will be freed immediately upon return (to UMF_kernel */
    /* and then to UMFPACK_numeric). */
}
