/* ========================================================================== */
/* === UMF_tuple_lengths ==================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* Determine the tuple list lengths, and the amount of memory required for */
/* them.  Return the amount of memory needed to store all the tuples. */
/* This routine assumes that the tuple lists themselves are either already */
/* deallocated, or will be shortly (so Row[ ].tlen and Col[ ].tlen are */
/* overwritten) */

#include "umf_internal.h"
#include "umf_build_tuples_usage.h"

GLOBAL Int UMF_tuple_lengths
(
    NumericType *Numeric,
    WorkType *Work,
    double *dusage
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int e, nrows, ncols, nel, i, *Rows, *Cols, row, col, n_row, n_col, *E,
	*Row_degree, *Row_tlen, *Col_degree, *Col_tlen, usage ;
    Element *ep ;
    Unit *p ;

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    E = Work->E ;
    Row_degree = Numeric->Rperm ;
    Col_degree = Numeric->Cperm ;
    Row_tlen   = Numeric->Uilen ;
    Col_tlen   = Numeric->Lilen ;
    n_row = Work->n_row ;
    n_col = Work->n_col ;
    nel = Work->nel ;

    DEBUG3 (("TUPLE_LENGTHS: n_row "ID" n_col "ID" nel "ID"\n",
	n_row, n_col, nel)) ;

    /* tuple list lengths already initialized to zero */

    /* ---------------------------------------------------------------------- */
    /* scan each element: count tuple list lengths (include element 0) */
    /* ---------------------------------------------------------------------- */

    for (e = 0 ; e <= nel ; e++)	/* for all elements, in any order */
    {
	if (E [e])
	{
#ifndef NDEBUG
	    UMF_dump_element (Numeric, Work, e, FALSE) ;
#endif
	    p = Numeric->Memory + E [e] ;
	    GET_ELEMENT_PATTERN (ep, p, Cols, Rows, ncols) ;
	    nrows = ep->nrows ;
	    ASSERT (e != 0) ;
	    for (i = 0 ; i < nrows ; i++)
	    {
		row = Rows [i] ;
		ASSERT (row == EMPTY || (row >= 0 && row < n_row)) ;
		if (row >= 0)
		{
		    ASSERT (NON_PIVOTAL_ROW (row)) ;
		    Row_tlen [row] ++ ;
		}
	    }
	    for (i = 0 ; i < ncols ; i++)
	    {
		col = Cols [i] ;
		ASSERT (col == EMPTY || (col >= 0 && col < n_col)) ;
		if (col >= 0)
		{
		    ASSERT (NON_PIVOTAL_COL (col)) ;
		    Col_tlen [col] ++ ;
		}
	    }
	}
    }

    /* note: tuple lengths are now modified, but the tuple lists are not */
    /* updated to reflect that fact. */

    /* ---------------------------------------------------------------------- */
    /* determine the required memory to hold all the tuple lists */
    /* ---------------------------------------------------------------------- */

    usage = UMF_build_tuples_usage (Col_tlen, Col_degree,
	Row_tlen, Row_degree, n_row, n_col, dusage) ;
    return (usage) ;
}

