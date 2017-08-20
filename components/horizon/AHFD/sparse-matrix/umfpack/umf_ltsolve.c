/* ========================================================================== */
/* === UMF_ltsolve ========================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*  Solves L'x = b or L.'x=b, where L is the lower triangular factor of a */
/*  matrix.  B is overwritten with the solution X. */
/*  Returns the floating point operation count */

#include "umf_internal.h"

GLOBAL double
#ifdef CONJUGATE_SOLVE
UMF_lhsolve			/* solve L'x=b  (complex conjugate transpose) */
#else
UMF_ltsolve			/* solve L.'x=b (array transpose) */
#endif
(
    NumericType *Numeric,
    Entry X [ ],		/* b on input, solution x on output */
    Int Pattern [ ]		/* a work array of size n */
)
{
    Int k, deg, *ip, j, row, *Lpos, *Lilen, kstart, kend, *Lip, llen,
	lp, pos, npiv ;
    Entry *xp, xk ;

    /* ---------------------------------------------------------------------- */

    if (Numeric->n_row != Numeric->n_col) return (0.) ;
    npiv = Numeric->npiv ;
    Lpos = Numeric->Lpos ;
    Lilen = Numeric->Lilen ;
    Lip = Numeric->Lip ;
    kstart = npiv ;

#ifndef NDEBUG
    DEBUG4 (("Ltsolve start:\n")) ;
    for (j = 0 ; j < Numeric->n_row ; j++)
    {
	DEBUG4 (("Ltsolve start "ID": ", j)) ;
	EDEBUG4 (X [j]) ;
	DEBUG4 (("\n")) ;
    }
#endif

    for (kend = npiv-1 ; kend >= 0 ; kend = kstart-1)
    {

	/* ------------------------------------------------------------------ */
	/* find the start of this Lchain */
	/* ------------------------------------------------------------------ */

	/* for (kstart = kend ; kstart >= 0 && Lip [kstart] > 0 ; kstart--) ; */
	kstart = kend ;
	while (kstart >= 0 && Lip [kstart] > 0)
	{
	    kstart-- ;
	}

	/* the Lchain goes from kstart to kend */

	/* ------------------------------------------------------------------ */
	/* scan the whole chain to find the pattern of the last column of L */
	/* ------------------------------------------------------------------ */

	deg = 0 ;
	DEBUG4 (("start of chain for column of L\n")) ;
	for (k = kstart ; k <= kend ; k++)
	{
	    ASSERT (k >= 0 && k < npiv) ;

	    /* -------------------------------------------------------------- */
	    /* make column k of L in Pattern [0..deg-1] */
	    /* -------------------------------------------------------------- */

	    /* remove pivot row */
	    pos = Lpos [k] ;
	    if (pos != EMPTY)
	    {
		DEBUG4 (("  k "ID" removing row "ID" at position "ID"\n",
		k, Pattern [pos], pos)) ;
		ASSERT (k != kstart) ;
		ASSERT (deg > 0) ;
		ASSERT (pos >= 0 && pos < deg) ;
		ASSERT (Pattern [pos] == k) ;
		Pattern [pos] = Pattern [--deg] ;
	    }

	    /* concatenate the pattern */
	    lp = Lip [k] ;
	    if (k == kstart)
	    {
		lp = -lp ;
	    }
	    ASSERT (lp > 0) ;
	    ip = (Int *) (Numeric->Memory + lp) ;
	    llen = Lilen [k] ;
	    for (j = 0 ; j < llen ; j++)
	    {
		row = *ip++ ;
		DEBUG4 (("  row "ID"  k "ID"\n", row, k)) ;
		ASSERT (row > k) ;
		Pattern [deg++] = row ;
	    }

	}
	/* Pattern [0..deg-1] is now the pattern of column kend */

	/* ------------------------------------------------------------------ */
	/* solve using this chain, in reverse order */
	/* ------------------------------------------------------------------ */

	DEBUG4 (("Unwinding Lchain\n")) ;
	for (k = kend ; k >= kstart ; k--)
	{

	    /* -------------------------------------------------------------- */
	    /* use column k of L */
	    /* -------------------------------------------------------------- */

	    ASSERT (k >= 0 && k < npiv) ;
	    lp = Lip [k] ;
	    if (k == kstart)
	    {
		lp = -lp ;
	    }
	    ASSERT (lp > 0) ;
	    llen = Lilen [k] ;
	    xp = (Entry *) (Numeric->Memory + lp + UNITS (Int, llen)) ;
	    xk = X [k] ;
	    for (j = 0 ; j < deg ; j++)
	    {
		DEBUG4 (("  row "ID"  k "ID" value", Pattern [j], k)) ;
		EDEBUG4 (*xp) ;
		DEBUG4 (("\n")) ;
#ifdef CONJUGATE_SOLVE
		/* xk -= X [Pattern [j]] * conjugate (*xp) ; */
		MULT_SUB_CONJ (xk, X [Pattern [j]], *xp) ;
#else
		/* xk -= X [Pattern [j]] * (*xp) ; */
		MULT_SUB (xk, X [Pattern [j]], *xp) ;
#endif
		xp++ ;
	    }
	    X [k] = xk ;

	    /* -------------------------------------------------------------- */
	    /* construct column k-1 of L */
	    /* -------------------------------------------------------------- */

	    /* un-concatenate the pattern */
	    deg -= llen ;

	    /* add pivot row */
	    pos = Lpos [k] ;
	    if (pos != EMPTY)
	    {
		DEBUG4 (("  k "ID" adding row "ID" at position "ID"\n",
		k, k, pos)) ;
		ASSERT (k != kstart) ;
		ASSERT (pos >= 0 && pos <= deg) ;
		Pattern [deg++] = Pattern [pos] ;
		Pattern [pos] = k ;
	    }
	}
    }

#ifndef NDEBUG
    for (j = 0 ; j < Numeric->n_row ; j++)
    {
	DEBUG4 (("Ltsolve done "ID": ", j)) ;
	EDEBUG4 (X [j]) ;
	DEBUG4 (("\n")) ;
    }
    DEBUG4 (("Ltsolve done.\n")) ;
#endif

    return (MULTSUB_FLOPS * ((double) Numeric->lnz)) ;
}

