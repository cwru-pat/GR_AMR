/* ========================================================================== */
/* === UMF_init_front ======================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

#include "umf_internal.h"

GLOBAL void UMF_init_front
(
    WorkType *Work
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int i, fsize, pivrow, pivcol, j, fnrows_max, fncols_max,
	row, col, *Frows, *Fcols, *Fcpos, *Frpos, fncols, fnrows, src, dest,
	*Wrow ;
    Entry *Fx, *F ;

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    /* current front is defined by pivot row and column */

    pivrow = Work->pivrow ;
    pivcol = Work->pivcol ;

    Frows = Work->Frows ;
    Fcols = Work->Fcols ;
    Frpos = Work->Frpos ;
    Fcpos = Work->Fcpos ;

    Work->fnpiv = 1 ;
    Work->fnzeros = 0 ;

    /* dynamic front dimensions, but fixed across one chain */
    fnrows_max = Work->fnrows_max ;
    fncols_max = Work->fncols_max ;

    fsize = fnrows_max * fncols_max ;

    /* ---------------------------------------------------------------------- */
    /* place pivot column pattern in frontal matrix */
    /* ---------------------------------------------------------------------- */

    if (Work->pivcol_in_front)
    {
	/* append the pivot column extension */
	/* note that all we need to do is increment the size, since the */
	/* candidate pivot column pattern is already in place in */
	/* Frows [0 ... Work->fnrows-1] (the old pattern), and */
	/* Frows [Work->fnrows ... Work->fnrows + Work->ccdeg - 1] (the new */
	/* pattern). */

	/* if both pivrow and pivcol are in front, then we extend the old one */
	/* in UMF_extend_front, rather than starting a new one here. */
	ASSERT (!Work->pivrow_in_front) ;

	dest = Work->fnrows  ;
	Work->fscan_row = dest ;	/* only scan the new rows */
	dest += Work->ccdeg ;
    }
    else
    {
	/* this is a completely new column */
	dest = 0 ;
	Work->fscan_row = 0 ;		/* scan all the rows */

	ASSERT (Work->Wcol == Work->Wm) ;
	for (i = 0 ; i < Work->ccdeg ; i++)
	{
	    row = Work->Wcol [i] ;
	    Frows [dest] = row ;
	    Frpos [row] = dest ;
	    dest++ ;
	}
    }

    Work->fnrows = dest ;

    /* place pivot row index into position */
    Frpos [pivrow] = fnrows_max - 1 ;

#ifndef NDEBUG
    DEBUG3 (("New Pivot col "ID" now in front, length "ID"\n",
	pivcol, Work->fnrows)) ;
    for (i = 0 ; i < Work->fnrows ; i++)
    {
	DEBUG4 (("row "ID" position "ID"\n", Frows [i], Frpos [Frows [i]])) ;
	ASSERT (Frpos [Frows [i]] == i) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* place pivot row pattern in frontal matrix */
    /* ---------------------------------------------------------------------- */

    if (Work->pivrow_in_front)
    {
	/* append the pivot row extension */
	src = Work->fncols ;
	dest = Work->fncols ;
	Work->fscan_col = dest ;	/* only scan the new columns */
    }
    else
    {
	/* this is a completely new row */
	src = 0 ;
	dest = 0 ;
	Work->fscan_col = 0 ;		/* scan all the columns */
    }

    Wrow = Work->Wrow ;
    for ( ; src < Work->rrdeg ; src++)
    {
	col = Wrow [src] ;
	if (col != pivcol)
	{
	    Fcols [dest] = col ;
	    Fcpos [col] = dest * fnrows_max ;
	    dest++ ;
	}
    }

    Work->fncols = dest ;

    /* place pivot column index into position */
    Fcpos [pivcol] = (fncols_max - 1) * fnrows_max ;

#ifndef NDEBUG
    DEBUG3 (("New Pivot row "ID" now in front, length "ID" fnrows_max "ID"\n",
		pivrow, Work->fncols, fnrows_max)) ;
    for (j = 0 ; j < Work->fncols ; j++)
    {
	DEBUG4 (("col "ID" position "ID"  ("ID")\n",
	    Fcols [j], Fcpos [Fcols [j]], j*fnrows_max)) ;
	ASSERT (Fcpos [Fcols [j]] == j * fnrows_max) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* clear the frontal matrix */
    /* ---------------------------------------------------------------------- */

    fncols = Work->fncols ;
    fnrows = Work->fnrows ;
    Fx = Work->Fx ;

    /* clear the contribution block and the pivot row */
    F = Fx ;
    for (j = 0 ; j < fncols ; j++)
    {
	for (i = 0 ; i < fnrows ; i++)
	{
	    ASSERT ((&F[i] >= Fx) && (&F[i] < Fx+fsize)) ;
	    CLEAR (F [i]) ;
	}
	ASSERT ((&F[fnrows_max-1] >= Fx) && (&F[fnrows_max-1] < Fx+fsize)) ;
	CLEAR (F [fnrows_max - 1]) ;
	F += fnrows_max ;
    }

    /* clear the pivot column (excl pivot itself) */
    F = Fx + Fcpos [pivcol] ;
    for (j = 0 ; j < fnrows ; j++)
    {
	ASSERT ((&F[j] >= Fx) && (&F[j] < Fx+fsize)) ;
	CLEAR (F [j]) ;
    }

    /* clear the pivot entry */
    j = fsize-1 ;
    DEBUG2 (("j "ID" fsize "ID"\n", j, fsize)) ;
    ASSERT ((&Fx[j] >= Fx) && (&Fx[j] < Fx+fsize)) ;
    CLEAR (Fx [j]) ;

    /* ---------------------------------------------------------------------- */
    /* current workspace usage: */
    /* ---------------------------------------------------------------------- */

    /* Fx [0..fnrows_max-1, 0..fncols_max-1]: */
    /*	space for the new frontal matrix. */
    /* Fx (i,j) is located at Fx [i+j*fnrows_max] */

    /* ---------------------------------------------------------------------- */
    /* make sure the pivot row and column are scanned in assembly phase */
    /* ---------------------------------------------------------------------- */

    Work->scan_pivrow = TRUE ;
    Work->scan_pivcol = TRUE ;

}
