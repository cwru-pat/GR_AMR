/* ========================================================================== */
/* === UMF_usolve =========================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*  solves Ux = b, where U is the upper triangular factor of a matrix. */
/*  B is overwritten with the solution X. */
/*  Returns the floating point operation count */

#include "umf_internal.h"

GLOBAL double UMF_usolve
(
    NumericType *Numeric,
    Entry X [ ],		/* b on input, solution x on output */
    Int Pattern [ ]		/* a work array of size n */
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int k, deg, j, *ip, col, *Upos, *Uilen, pos,
	*Uip, n, ulen, up, newUchain, npiv ;
    Entry *xp, xk, *D ;
    double flops ;

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    if (Numeric->n_row != Numeric->n_col) return (0.) ;
    n = Numeric->n_row ;
    npiv = Numeric->npiv ;
    Upos = Numeric->Upos ;
    Uilen = Numeric->Uilen ;
    Uip = Numeric->Uip ;
    D = Numeric->D ;

    /* start by counting the division by D [0..n-1] */
    flops = n ;

#ifndef NDEBUG
    DEBUG4 (("Usolve start:  npiv = "ID" n = "ID"\n", npiv, n)) ;
    for (j = 0 ; j < n ; j++)
    {
	DEBUG4 (("Usolve start "ID": ", j)) ;
	EDEBUG4 (X [j]) ;
	DEBUG4 (("\n")) ;
    }
#endif

    /* handle the singular part of D, up to just before the last pivot */
    flops += (n-npiv) ;
    for (k = n-1 ; k >= npiv ; k--)
    {
	/* Note that D [k] is identically zero. */
	ASSERT (IS_ZERO (D [k])) ;
	xk = X [k] ;
	/* X [k] = xk / D [k] ; */
	DIV (X [k], xk, D [k]) ;
    }

    deg = Numeric->ulen ;
    if (deg > 0)
    {
	/* make last pivot row of U (singular matrices only) */
	for (j = 0 ; j < deg ; j++)
	{
	    DEBUG1 (("Last row of U: j="ID"\n", j)) ;
	    DEBUG1 (("Last row of U: Upattern[j]="ID"\n",
		Numeric->Upattern [j]) );
	    Pattern [j] = Numeric->Upattern [j] ;
	}
    }

    for (k = npiv-1 ; k >= 0 ; k--)
    {

	/* ------------------------------------------------------------------ */
	/* use row k of U */
	/* ------------------------------------------------------------------ */

	up = Uip [k] ;
	ulen = Uilen [k] ;
	newUchain = (up < 0) ;
	if (newUchain)
	{
	    up = -up ;
	    xp = (Entry *) (Numeric->Memory + up + UNITS (Int, ulen)) ;
	}
	else
	{
	    xp = (Entry *) (Numeric->Memory + up) ;
	}

	xk = X [k] ;
	flops += 2*deg ;
	for (j = 0 ; j < deg ; j++)
	{
	    DEBUG4 (("  k "ID" col "ID" value", k, Pattern [j])) ;
	    EDEBUG4 (*xp) ;
	    DEBUG4 (("\n")) ;
	    /* xk -= X [Pattern [j]] * (*xp) ; */
	    MULT_SUB (xk, X [Pattern [j]], *xp) ;
	    xp++ ;
	}

	/* Go ahead and divide by zero if D [k] is zero */
	/* X [k] = xk / D [k] ; */
	DIV (X [k], xk, D [k]) ;

	/* ------------------------------------------------------------------ */
	/* make row k-1 of U in Pattern [0..deg-1] */
	/* ------------------------------------------------------------------ */

	if (newUchain)
	{
	    /* next row is a new Uchain */
	    deg = ulen ;
	    ASSERT (IMPLIES (k == 0, deg == 0)) ;
	    DEBUG4 (("end of chain for row of U "ID" deg "ID"\n", k-1, deg)) ;
	    ip = (Int *) (Numeric->Memory + up) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		col = *ip++ ;
		DEBUG4 (("  k "ID" col "ID"\n", k-1, col)) ;
		ASSERT (k <= col) ;
		Pattern [j] = col ;
	    }
	}
	else
	{
	    deg -= ulen ;
	    DEBUG4 (("middle of chain for row of U "ID" deg "ID"\n", k, deg)) ;
	    ASSERT (deg >= 0) ;
	    pos = Upos [k] ;
	    if (pos != EMPTY)
	    {
		/* add the pivot column */
		DEBUG4 (("k "ID" add pivot entry at pos "ID"\n", k, pos)) ;
		ASSERT (pos >= 0 && pos <= deg) ;
		Pattern [deg++] = Pattern [pos] ;
		Pattern [pos] = k ;
	    }
	}
    }

#ifndef NDEBUG
    for (j = 0 ; j < n ; j++)
    {
	DEBUG4 (("Usolve done "ID": ", j)) ;
	EDEBUG4 (X [j]) ;
	DEBUG4 (("\n")) ;
    }
    DEBUG4 (("Usolve done.\n")) ;
#endif

    return (DIV_FLOPS * ((double) n) + MULTSUB_FLOPS * ((double) Numeric->unz));
}

