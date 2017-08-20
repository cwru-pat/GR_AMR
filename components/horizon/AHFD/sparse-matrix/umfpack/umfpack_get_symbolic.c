/* ========================================================================== */
/* === UMFPACK_get_symbolic ================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Gets the symbolic information held in the Symbolic object.
    See umfpack_get_symbolic.h for a more detailed description.
*/

#include "umf_internal.h"
#include "umf_valid_symbolic.h"

/* ========================================================================== */

GLOBAL Int UMFPACK_get_symbolic
(
    Int *p_n_row,
    Int *p_n_col,
    Int *p_nz,
    Int *p_nfr,
    Int *p_nchains,
    Int Ptree [ ],
    Int Qtree [ ],
    Int Front_npivcol [ ],
    Int Front_parent [ ],
    Int Front_1strow [ ],
    Int Front_leftmostdesc [ ],
    Int Chain_start [ ],
    Int Chain_maxrows [ ],
    Int Chain_maxcols [ ],
    void *SymbolicHandle
)
{
    SymbolicType *Symbolic ;
    Int k, n_row, n_col, nfr, nchains, *p ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    Symbolic = (SymbolicType *) SymbolicHandle ;
    if (!UMF_valid_symbolic (Symbolic))
    {
	return (UMFPACK_ERROR_invalid_Symbolic_object) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get contents of Symbolic */
    /* ---------------------------------------------------------------------- */

    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;
    nfr = Symbolic->nfr ;
    nchains = Symbolic->nchains ;

    if (p_n_row)
    {
	*p_n_row = n_row ;
    }

    if (p_n_col)
    {
	*p_n_col = n_col ;
    }

    if (p_nz)
    {
	*p_nz = Symbolic->nz ;
    }

    if (p_nfr)
    {
	*p_nfr = nfr ;
    }

    if (p_nchains)
    {
	*p_nchains = nchains ;
    }

    if (Ptree)
    {
	p = Symbolic->Rperm_init ;
	for (k = 0 ; k < n_row ; k++)
	{
	    Ptree [k] = p [k] ;
	}
    }

    if (Qtree)
    {
	p = Symbolic->Cperm_init ;
	for (k = 0 ; k < n_col ; k++)
	{
	    Qtree [k] = p [k] ;
	}
    }

    if (Front_npivcol)
    {
	p = Symbolic->Front_npivcol ;
	for (k = 0 ; k <= nfr ; k++)
	{
	    Front_npivcol [k] = p [k] ;
	}
    }

    if (Front_parent)
    {
	p = Symbolic->Front_parent ;
	for (k = 0 ; k <= nfr ; k++)
	{
	    Front_parent [k] = p [k] ;
	}
    }

    if (Front_1strow)
    {
	p = Symbolic->Front_1strow ;
	for (k = 0 ; k <= nfr ; k++)
	{
	    Front_1strow [k] = p [k] ;
	}
    }

    if (Front_leftmostdesc)
    {
	p = Symbolic->Front_leftmostdesc ;
	for (k = 0 ; k <= nfr ; k++)
	{
	    Front_leftmostdesc [k] = p [k] ;
	}
    }

    if (Chain_start)
    {
	p = Symbolic->Chain_start ;
	for (k = 0 ; k <= nchains ; k++)
	{
	    Chain_start [k] = p [k] ;
	}
    }

    if (Chain_maxrows)
    {
	p = Symbolic->Chain_maxrows ;
	for (k = 0 ; k < nchains ; k++)
	{
	    Chain_maxrows [k] = p [k] ;
	}
	Chain_maxrows [nchains] = 0 ;
    }

    if (Chain_maxcols)
    {
	p = Symbolic->Chain_maxcols ;
	for (k = 0 ; k < nchains ; k++)
	{
	    Chain_maxcols [k] = p [k] ;
	}
	Chain_maxcols [nchains] = 0 ;
    }

    return (UMFPACK_OK) ;
}
