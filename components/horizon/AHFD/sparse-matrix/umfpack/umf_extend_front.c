/* ========================================================================== */
/* === UMF_extend_front ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* Called by kernel. */

#include "umf_internal.h"

GLOBAL void UMF_extend_front
(
    WorkType *Work
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Entry *Fx, *Flast, *Fd, *Fs, *Fcol, *Frow, temp, *Flrow, *Fu ;
    Int j, i, fnpiv, *Frows, row, col, i2, *Wrow, *Wcol,
	*Frpos, *Fcpos, *Fcols, pivcol, pivrow, fnrows_extended, rrdeg, ccdeg,
	fncols_extended, fnrows_max, fncols_max, fnrows, fncols,
	fspos, fdpos, fs, j2, j3, col2, row2, fncols_orig ;
#ifndef NDEBUG
    Int debug_ok ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    Frows = Work->Frows ;
    Frpos = Work->Frpos ;
    Fcols = Work->Fcols ;
    Fcpos = Work->Fcpos ;
    Fx = Work->Fx ;
    pivcol = Work->pivcol ;
    pivrow = Work->pivrow ;
    fnrows_max = Work->fnrows_max ;
    fncols_max = Work->fncols_max ;
    fnrows = Work->fnrows ;
    fncols = Work->fncols ;
    fncols_orig = fncols ;
    rrdeg = Work->rrdeg ;
    ccdeg = Work->ccdeg ;

#ifndef NDEBUG
    DEBUG2 (("EXTEND FRONT\n")) ;
    if (fncols == 0 || fnrows == 0)
    {
	DEBUG2 (("Extending empty front "ID" "ID"\n", fnrows,fncols)) ;
    }
    DEBUG6 (("Pivot row, before shift and extension: "ID"\n", fncols)) ;
    for (j = 0 ; j < fncols ; j++)
    {
	DEBUG7 ((" "ID" "ID" "ID"  "ID"\n",
	    j, Fcols [j], Fcpos [Fcols [j]], j < fncols)) ;
	ASSERT (Fcpos [Fcols [j]] == j * fnrows_max) ;
    }
    DEBUG6 (("Pivot col, before shift and extension: "ID"\n", fnrows)) ;
    for (i = 0 ; i < fnrows ; i++)
    {
	DEBUG7 ((" "ID" "ID" "ID"  "ID"\n",
	    i, Frows [i], Frpos [Frows [i]], i < fnrows)) ;
	ASSERT (Frpos [Frows [i]] == i) ;
    }
    if (Work->pivcol_in_front)
    {
	DEBUG6 (("Extended part of pivot col, before shift/extension: "ID"\n",
	    fnrows)) ;
	for (i = fnrows ; i < fnrows + ccdeg ; i++)
	{
	    DEBUG7 ((" "ID" "ID" "ID"  "ID"\n",
		i, Frows [i], Frpos [Frows [i]], i < fnrows)) ;
	    ASSERT (Frpos [Frows [i]] == i) ;
	}
    }
    DEBUG2 (("Work->fnpiv "ID"\n", Work->fnpiv)) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* pivot column already updated */
    /* ---------------------------------------------------------------------- */

    fnpiv = Work->fnpiv ;

    /* ====================================================================== */
    /* === Shift pivot row and column ======================================= */
    /* ====================================================================== */

    /* ---------------------------------------------------------------------- */
    /* move the pivot column into place */
    /* ---------------------------------------------------------------------- */

    fdpos = (fncols_max - fnpiv - 1) * fnrows_max ;
    Fd = Fx + fdpos ;	/* Fd: destination of pivot column */

#ifndef NDEBUG
    /*
    DEBUG7 (("Complete frontal matrix prior to pivcol swap (incl unused):\n")) ;
    UMF_dump_dense (Fx, fnrows_max, fnrows_max, fncols_max) ;
    */
#endif

    if (Work->pivcol_in_front)
    {

	fspos = Fcpos [pivcol] ;
	fs = fspos / fnrows_max ;

	/* Fs: current position of pivot column */
	Fs = Fx + fspos ;

	/* Flast: position of last column in front */
	Flast = Fx + (fncols - 1) * fnrows_max ;

	/* ------------------------------------------------------------------ */
	/* pivot column is in current front - shift into place */
	/* ------------------------------------------------------------------ */

	DEBUG6 (("Swap/shift pivot column in front\n")) ;
	DEBUG6 (("fspos: "ID" flpos: "ID" fdpos: "ID"\n",
	    fspos, (fncols - 1) * fnrows_max, fdpos)) ;

	if (Flast != Fd)
	{
	    if (Fs == Flast)
	    {
		/* ---------------------------------------------------------- */
		/* move Fs => Fd */
		/* ---------------------------------------------------------- */

		DEBUG6 (("col case 1\n")) ;

		/* column of the contribution block: */
		for (i = 0 ; i < fnrows ; i++)
		{
		    Fd [i] = Fs [i] ;
		}

#ifndef NDEBUG
		/* column of the U2 block */
		for (i = fnrows_max - fnpiv ; i < fnrows_max ; i++)
		{
		    ASSERT (IS_ZERO (Fs [i])) ;
		}
#endif

	    }
	    else
	    {
		/* ---------------------------------------------------------- */
		/* move Fs => Fd */
		/* move Flast => Fs */
		/* ---------------------------------------------------------- */

		DEBUG6 (("col case 2\n")) ;

		/* column of the contribution block: */
		for (i = 0 ; i < fnrows ; i++)
		{
		    Fd [i] = Fs [i] ;
		    Fs [i] = Flast [i] ;
		}
		/* column of the U2 block */
		for (i = fnrows_max - fnpiv ; i < fnrows_max ; i++)
		{
		    ASSERT (IS_ZERO (Fs [i])) ;
		    Fs [i] = Flast [i] ;
		}
	    }
	}
	else if (Fs != Fd)
	{
	    /* -------------------------------------------------------------- */
	    /* swap Fs <=> Fd */
	    /* -------------------------------------------------------------- */

	    DEBUG6 (("col case 3\n")) ;

	    /* column of the contribution block: */
	    for (i = 0 ; i < fnrows ; i++)
	    {
		temp = Fd [i] ;
		Fd [i] = Fs [i] ;
		Fs [i] = temp ;
	    }
	    /* column of the U2 block */
	    for (i = fnrows_max - fnpiv ; i < fnrows_max ; i++)
	    {
		ASSERT (IS_ZERO (Fs [i])) ;
		Fs [i] = Fd [i] ;
	    }
	}

	/* move column Flast to Fs in the Fcols pattern */
	col2 = Fcols [fncols - 1] ;
	Fcols [fs] = col2 ;
	Fcpos [col2] = fspos ;

	/* one less column in the contribution block */
	fncols = --(Work->fncols) ;

    }
    else
    {
	/* ------------------------------------------------------------------ */
	/* pivot column is not in front - tack onto L block */
	/* ------------------------------------------------------------------ */

	DEBUG6 (("Pivot column not in front\ncol case 5\n")) ;

	/* column of L */
	for (i = 0 ; i < fnrows ; i++)
	{
	    CLEAR (Fd [i]) ;
	}
	/* column of U2 */
	for (i = fnrows_max - fnpiv ; i < fnrows_max ; i++)
	{
	    CLEAR (Fd [i]) ;
	}

    }

    /* move pivot column to Fd */
    Fcpos [pivcol] = fdpos ;

    /* scan starts at the first new column in Fcols */
    /* also scan the pivot column if it was not in the front */
    Work->fscan_col = fncols ;

#ifndef NDEBUG
    DEBUG6 (("Pivot row, after shift but before extension: "ID"\n", fncols)) ;
    for (j = 0 ; j < fncols ; j++)
    {
	DEBUG7 ((" "ID" "ID" "ID"  "ID"\n",
	    j, Fcols [j], Fcpos [Fcols [j]], j < fncols)) ;
	ASSERT (Fcpos [Fcols [j]] == j * fnrows_max) ;
    }
    /*
    DEBUG7 (("Complete frontal matrix after pivot col swap (incl unused):\n")) ;
    UMF_dump_dense (Fx, fnrows_max, fnrows_max, fncols_max) ;
    */
#endif

    /* ---------------------------------------------------------------------- */
    /* move pivot row into place */
    /* ---------------------------------------------------------------------- */

    fdpos = fnrows_max - fnpiv - 1 ;
    Fd = Fx + fdpos ;	/* Fd: destination of pivot row */

    if (Work->pivrow_in_front)
    {

	fspos = Frpos [pivrow] ;

	/* Fs: current position of pivot column in front */
	Fs = Fx + fspos ;

	/* Flast: position of last row in front */
	Flast = Fx + (fnrows - 1) ;

	/* ------------------------------------------------------------------ */
	/* pivot row is in current front - shift into place */
	/* ------------------------------------------------------------------ */

	DEBUG6 (("Swap/shift pivot row in front:\n")) ;
	DEBUG6 (("fspos: "ID" flpos: "ID" fdpos: "ID"\n",
	    fspos, fnrows-1, fdpos)) ;

	if (Flast != Fd)
	{
	    if (Fs == Flast)
	    {
		/* ---------------------------------------------------------- */
		/* move Fs => Fd */
		/* ---------------------------------------------------------- */

		DEBUG6 (("row case 1\n")) ;

		/* row of the contribution block: */
		j2 = fncols * fnrows_max ;
		for (j = 0 ; j < j2 ; j += fnrows_max)
		{
		    Fd [j] = Fs [j] ;
		}

		/* row of the L2 block: */
		j2 = fncols_max	* fnrows_max ;
		j = (fncols_max - fnpiv	- 1) * fnrows_max ;
		for ( ; j < j2 ; j += fnrows_max)
		{
		    Fd [j] = Fs [j] ;
		}

	    }
	    else
	    {

		/* ---------------------------------------------------------- */
		/* move Fs => Fd */
		/* move Flast => Fs */
		/* ---------------------------------------------------------- */

		DEBUG6 (("row case 2\n")) ;

		/* row of the contribution block: */
		j2 = fncols * fnrows_max ;
		for (j = 0 ; j < j2 ; j += fnrows_max)
		{
		    Fd [j] = Fs [j] ;
		    Fs [j] = Flast [j] ;
		}
		/* row of the L2 block: */
		j2 = fncols_max * fnrows_max ;
		j = (fncols_max - fnpiv - 1) * fnrows_max ;
		for ( ; j < j2 ; j += fnrows_max)
		{
		    Fd [j] = Fs [j] ;
		    Fs [j] = Flast [j] ;
		}
	    }
	}
	else if (Fs != Fd)
	{

	    /* -------------------------------------------------------------- */
	    /* swap Fs <=> Fd */
	    /* -------------------------------------------------------------- */

	    DEBUG6 (("row case 3\n")) ;

	    /* row of the contribution block: */
	    j2 = fncols * fnrows_max ;
	    for (j = 0 ; j < j2 ; j += fnrows_max)
	    {
		temp = Fd [j] ;
		Fd [j] = Fs [j] ;
		Fs [j] = temp ;
	    }
	    /* row of the L2 block: */
	    j2 = fncols_max * fnrows_max ;
	    j = (fncols_max - fnpiv - 1) * fnrows_max ;
	    for ( ; j < j2 ; j += fnrows_max)
	    {
		temp = Fd [j]  ;
		Fd [j] = Fs [j] ;
		Fs [j] = temp ;
	    }
	}

	/* move row Flast to Fs in the Frows pattern */
	row2 = Frows [fnrows-1] ;
	Frows [fspos] = row2 ;
	Frpos [row2] = fspos ;


	if (Work->pivcol_in_front && ccdeg > 0)
	{

	    /* move row Fe to Flast in the extended Frows pattern */
	    row2 = Frows [fnrows + ccdeg - 1] ;
	    Frows [fnrows-1] = row2 ;
	    Frpos [row2] = fnrows-1 ;
	}

	/* one less row in the contribution block */
	fnrows = --(Work->fnrows) ;

	/* ------------------------------------------------------------------ */
	/* update pivot row */
	/* ------------------------------------------------------------------ */

	if (fnpiv > 0 && fncols > 0)
	{
	    DEBUG6 (("Update pivot row (but not pivot entry itself):\n")) ;
	    Fu = Fd + 1 ;
	    Flrow = Fd + (fncols_max - fnpiv) * fnrows_max ;

#ifdef USE_NO_BLAS

	    /* no BLAS available - use plain C code instead */
	    j2 = 0 ;
	    for (j = 0 ; j < fncols ; j++)
	    {
		i2 = 0 ;
		for (i = 0 ; i < fnpiv ; i++)
		{
		    /* Fd [j2] -= Fu [i+j*fnrows_max] * Flrow [i2] ; */
		    MULT_SUB (Fd [j2], Fu [i+j*fnrows_max], Flrow [i2]) ;
		    i2 += fnrows_max ;
		}
		j2 += fnrows_max ;
	    }

#else

	    BLAS_GEMV_ROW (fnpiv, fncols, Fu, Flrow, Fd, fnrows_max) ;

#endif	/* USE_NO_BLAS */

	}

    }
    else
    {
	/* ------------------------------------------------------------------ */
	/* pivot row is not in front - tack onto U block */
	/* ------------------------------------------------------------------ */

	DEBUG6 (("Pivot row not in current front\nrow case 5\n")) ;

	/* row of U */
	j2 = fncols * fnrows_max ;
	for (j = 0 ; j < j2 ; j += fnrows_max)
	{
	    CLEAR (Fd [j]) ;
	}
	/* row of L2 */
	j2 = fncols_max * fnrows_max ;
	j = (fncols_max - fnpiv - 1) * fnrows_max ;
	for ( ; j < j2 ; j += fnrows_max)
	{
	    CLEAR (Fd [j]) ;
	}

    }

    /* move pivot row to Fd */
    Frpos [pivrow] = fdpos ;

    /* scan1 starts at the first new row in Frows */
    /* also scan the pivot row if it was not in the front */
    Work->fscan_row = fnrows ;

#ifndef NDEBUG
    debug_ok = TRUE ;
    DEBUG6 (("Pivot col, before shift but before extension: "ID"\n", fnrows)) ;
    for (i = 0 ; i < fnrows ; i++)
    {
	DEBUG7 ((" "ID" "ID" "ID"  "ID"\n",
	    i, Frows [i], Frpos [Frows [i]], i < fnrows)) ;
	debug_ok = debug_ok && (Frpos [Frows [i]] == i) ;
    }
    ASSERT (debug_ok) ;
#endif

    /* ====================================================================== */
    /* === EXTEND PATTERN OF FRONT ========================================== */
    /* ====================================================================== */

    /* ---------------------------------------------------------------------- */
    /* extend row pattern of the front with the new pivot column extension */
    /* ---------------------------------------------------------------------- */

    fnrows_extended = fnrows ;
    fncols_extended = fncols ;

    if (Work->pivcol_in_front)
    {
	/* extended pattern and position already in Frows and Frpos */
	fnrows_extended += ccdeg ;

#ifndef NDEBUG
    debug_ok = TRUE ;
    for (i = fnrows ; i < fnrows + ccdeg ; i++)
    {
	row = Frows [i] ;
	DEBUG2 ((" row:: "ID" (ext)\n", row)) ;
	debug_ok = debug_ok && (Frpos [row] == i) && (row != pivrow) ;
    }
    ASSERT (debug_ok) ;
#endif

    }
    else
    {
	/* extended pattern is in Wcol, not yet in the front */
	Wcol = Work->Wcol ;
	ASSERT (Wcol == Work->Wm) ;
	for (i = 0 ; i < ccdeg ; i++)
	{
	    row = Wcol [i] ;
	    DEBUG2 ((" row:: "ID" (ext)\n", row)) ;
	    ASSERT (Frpos [row] == EMPTY);
	    ASSERT (row != pivrow) ;
	    Frows [fnrows_extended] = row ;
	    Frpos [row] = fnrows_extended ;
	    fnrows_extended++ ;
	}
    }

    ASSERT (fnpiv + fnrows_extended <= fnrows_max) ;

    /* ---------------------------------------------------------------------- */
    /* extend the column pattern of the front with the new pivot row */
    /* ---------------------------------------------------------------------- */

    if (Work->pivrow_in_front)
    {
	if (Work->pivcol_in_front)
	{

	    /* fill in the hole made when the pivot column was removed */
	    if (rrdeg > fncols_orig)
	    {
		Fcols [fncols] = Fcols [--rrdeg] ;
		fncols_extended = rrdeg ;
		for (i = fncols ; i < rrdeg ; i++)
		{
#ifndef NDEBUG
		    col = Fcols [i] ;
		    ASSERT (col != pivcol) ;
		    DEBUG2 ((" col:: "ID" (ext)\n", col)) ;
		    ASSERT (Fcpos [col] < 0) ;
#endif
		    Fcpos [Fcols [i]] = i * fnrows_max ;
		}
	    }
	}
	else
	{
	    Wrow = Work->Wrow ;
	    for (i = fncols_orig ; i < rrdeg ; i++)
	    {
		col = Wrow [i] ;
		if (col != pivcol)
		{
		    DEBUG2 ((" col:: "ID" (ext)\n", col)) ;
		    ASSERT (Fcpos [col] < 0) ;
		    Fcols [fncols_extended] = col ;
		    Fcpos [col] = fncols_extended * fnrows_max ;
		    fncols_extended++ ;
		}
	    }
	}
    }
    else
    {
	Wrow = Work->Wrow ;
	for (i = 0 ; i < rrdeg ; i++)
	{
	    col = Wrow [i] ;
	    if (col != pivcol && Fcpos [col] < 0)
	    {
		DEBUG2 ((" col:: "ID" (ext)\n", col)) ;
		Fcols [fncols_extended] = col ;
		Fcpos [col] = fncols_extended * fnrows_max ;
		fncols_extended++ ;
	    }
	}
    }


#ifndef NDEBUG
    ASSERT (fnpiv + fncols_extended <= fncols_max) ;
    DEBUG6 (("Pivot row, after shift and extension: "ID" "ID"\n",
    fncols,fncols_extended)) ;
    for (j = 0 ; j < fncols_extended ; j++)
    {
	DEBUG7 ((" "ID" "ID" "ID"  "ID"\n",
	    j, Fcols [j], Fcpos [Fcols [j]], j < fncols)) ;
	ASSERT (Fcpos [Fcols [j]] == j * fnrows_max) ;
    }
    DEBUG6 (("Pivot col, after shift and extension: "ID" "ID"\n",
    fnrows,fnrows_extended)) ;
    for (i = 0 ; i < fnrows_extended ; i++)
    {
	DEBUG7 ((" "ID" "ID" "ID"  "ID"\n",
	    i, Frows [i], Frpos [Frows [i]], i < fnrows)) ;
	ASSERT (Frpos [Frows [i]] == i) ;
    }
    /*
    DEBUG7 (("Complete frontal matrix after all swaps (incl unused):\n")) ;
    UMF_dump_dense (Fx, fnrows_max, fnrows_max, fncols_max) ;
    */
#endif

    /* ====================================================================== */
    /* Prepare for degree update and next local pivot search */
    /* ====================================================================== */

#ifndef NDEBUG
    DEBUG6 (("JUST BEFORE SCAN3A/4A:\nPivot row pattern:\n")) ;
    for (j = 0 ; j < fncols_extended ; j++)
    {
	DEBUG7 ((ID" "ID" "ID" "ID"\n",
	    j, Fcols [j], Fcpos [Fcols [j]], j < fncols)) ;
	ASSERT (Fcpos [Fcols [j]] == j * fnrows_max) ;
    }
    DEBUG6 (("Pivot col pattern:\n")) ;
    for (i = 0 ; i < fnrows_extended ; i++)
    {
	DEBUG7 ((ID" "ID" "ID" "ID"\n",
	    i, Frows [i], Frpos [Frows [i]], i < fnrows)) ;
	ASSERT (Frpos [Frows [i]] == i) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* Finished with step fnpiv (except for assembly and scale of pivot col) */
    /* ---------------------------------------------------------------------- */

    fnpiv = ++(Work->fnpiv) ;

#ifndef NDEBUG
    DEBUG6 (("EXT: pivot row pattern:  len="ID"\n", fncols_extended)) ;
    for (j = 0 ; j < fncols_extended ; j++) DEBUG7 ((ID"\n", Fcols [j])) ;
    DEBUG6 (("EXT: pivot col pattern:  len="ID"\n", fnrows_extended)) ;
    for (j = 0 ; j < fnrows_extended ; j++) DEBUG7 ((ID"\n", Frows [j])) ;
#endif

    /* ====================================================================== */
    /* === EXTEND NUMERICAL FRONT =========================================== */
    /* ====================================================================== */

    /* ---------------------------------------------------------------------- */
    /* Zero the newly extended frontal matrix */
    /* ---------------------------------------------------------------------- */

    Fcol = Fx + fncols * fnrows_max ;
    i2 = fnrows_max - fnpiv ;
    for (j = fncols ; j < fncols_extended ; j++)
    {
	/* zero the new columns in the contribution block: */
	for (i = 0 ; i < fnrows_extended ; i++)
	{
	    CLEAR (Fcol [i]) ;
	}
	/* zero the new columns in U block: */
	for (i = i2 ; i < fnrows_max ; i++)
	{
	    CLEAR (Fcol [i]) ;
	}
	Fcol += fnrows_max ;
    }

    Frow = Fx + fnrows ;
    j3 = fncols_max - fnpiv ;
    for (i = fnrows ; i < fnrows_extended ; i++)
    {
	/* zero the new rows in the contribution block: */
	for (j = 0 ; j < fncols ; j++)
	{
	    CLEAR (Frow [j * fnrows_max]) ;
	}
	/* zero the new rows in L block: */
	for (j = j3 ; j < fncols_max ; j++)
	{
	    CLEAR (Frow [j * fnrows_max]) ;
	}
	Frow++ ;
    }

    /* ---------------------------------------------------------------------- */
    /* finalize extended row and column pattern of the frontal matrix */
    /* ---------------------------------------------------------------------- */

    Work->fnrows = fnrows_extended ;
    Work->fncols = fncols_extended ;

    Work->scan_pivcol = !Work->pivcol_in_front ;
    Work->scan_pivrow = !Work->pivrow_in_front ;

}
