/* ========================================================================== */
/* === UMF_kernel =========================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Primary factorization routine.   Called by UMFPACK_numeric.
    Returns:
	UMFPACK_OK if successful,
	UMFPACK_ERROR_out_of_memory if out of memory, or
	UMFPACK_ERROR_different_pattern if pattern of matrix (Ap and/or Ai)
	   has changed since the call to UMFPACK_*symbolic.
*/

#include "umf_internal.h"
#include "umf_init_front.h"
#include "umf_assemble.h"
#include "umf_scale_column.h"
#include "umf_local_search.h"
#include "umf_create_element.h"
#include "umf_extend_front.h"
#include "umf_blas3_update.h"
#include "umf_kernel_init.h"
#include "umf_kernel_wrapup.h"


GLOBAL Int UMF_kernel
(
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int j, frontid1, frontid2, chain, nchains, *Chain_start, status,
	*Chain_maxrows, *Chain_maxcols, *Front_npivcol, jmax, nb ;

    /* ---------------------------------------------------------------------- */
    /* initialize memory space and load the matrix */
    /* ---------------------------------------------------------------------- */

    if (!UMF_kernel_init (Ap, Ai, Ax,
#ifdef COMPLEX
	Az,
#endif
	Numeric, Work, Symbolic))
    {
	/* UMF_kernel_init is guaranteed to succeed, since UMFPACK_numeric */
	/* either allocates enough space or if not, UMF_kernel does not get */
	/* called.  So running out of memory here is a fatal error, and means */
	/* that the user changed Ap and/or Ai since the call to */
	/* UMFPACK_*symbolic. */
	return (UMFPACK_ERROR_different_pattern) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get the symbolic factorization */
    /* ---------------------------------------------------------------------- */

    nchains = Symbolic->nchains ;
    Chain_start = Symbolic->Chain_start ;
    Chain_maxrows = Symbolic->Chain_maxrows ;
    Chain_maxcols = Symbolic->Chain_maxcols ;
    Front_npivcol = Symbolic->Front_npivcol ;
    DEBUG0 (("Starting Kernel, nchains "ID"\n", nchains)) ;
    nb = Symbolic->nb ;
    Work->nextcand = 0 ;

    /* ---------------------------------------------------------------------- */
    /* factorize each chain of frontal matrices */
    /* ---------------------------------------------------------------------- */

    for (chain = 0 ; chain < nchains ; chain++)
    {
	ASSERT (Work->fnrows == 0 && Work->fncols == 0 && Work->fnpiv == 0) ;
	frontid1 = Chain_start [chain] ;
	frontid2 = Chain_start [chain+1] - 1 ;
	Work->fnrows_max = Chain_maxrows [chain] ;
	Work->fncols_max = Chain_maxcols [chain] ;
	DEBUG2 (("Starting chain "ID". start "ID" end "ID" maxrows "ID
	    " maxcols "ID"\n", chain, frontid1, frontid2, Work->fnrows_max,
	    Work->fncols_max)) ;

	/* ------------------------------------------------------------------ */
	/* factorize each front in the chain */
	/* ------------------------------------------------------------------ */

	for (Work->frontid = frontid1 ; Work->frontid <= frontid2 ; Work->frontid++)
	{

	    /* -------------------------------------------------------------- */
	    /* there are no 0-by-c or r-by-0 fronts, where c and r are > 0 */
	    /* -------------------------------------------------------------- */

	    /* a front is either 0-by-0, or r-by-c */
	    DEBUG2 (("\n\n::: "ID" : Npiv: "ID" size "ID"-by-"ID"\n",
	    	Work->frontid, Work->npiv,
		Work->fnrows, Work->fncols)) ;
	    ASSERT ((Work->fnrows == 0 && Work->fncols == 0)
		  ||(Work->fnrows != 0 && Work->fncols != 0)) ;

	    /* -------------------------------------------------------------- */
	    /* Initialize the pivot column candidate set */
	    /* -------------------------------------------------------------- */

	    /* next pivot (in range 0 to n-1) */
	    Work->ncand = Front_npivcol [Work->frontid] ;
	    jmax = MIN (MAX_CANDIDATES, Work->ncand) ;
	    for (j = 0 ; j < jmax ; j++)
	    {
		Work->Candidates [j] = Work->nextcand++ ;
		DEBUG3 ((""ID" Candidate: "ID"\n", j, Work->Candidates [j])) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* Assemble and factorize the current frontal matrix */
	    /* -------------------------------------------------------------- */

	    while (Work->ncand > 0)
	    {

		/* ---------------------------------------------------------- */
		/* get the pivot row and column */
		/* ---------------------------------------------------------- */

		status = UMF_local_search (Numeric, Work, Symbolic) ;
		if (status == UMFPACK_ERROR_different_pattern)
		{
		    return (UMFPACK_ERROR_different_pattern) ;
		}
		if (status == UMFPACK_WARNING_singular_matrix)
		{
		    /* no pivot found, try again */
		    DEBUG0 (("No pivot found; try again, ncand: "ID"\n",
		    	Work->ncand)) ;
		    continue ;
		}

		/* ---------------------------------------------------------- */
		/* extend the frontal matrix, or start a new one */
		/* ---------------------------------------------------------- */

		if (Work->do_extend)
		{
		    /* apply pending updates if front would grow too much */
		    if (Work->do_update)
		    {
			UMF_blas3_update (Work) ;
		    }
		    /* extend the current front */
		    UMF_extend_front (Work) ;
		}
		else
		{
		    /* finish the current front (if any) and start a new one */
		    if (!UMF_create_element (Numeric, Work))
		    {
			return (UMFPACK_ERROR_out_of_memory) ;
		    }
		    UMF_init_front (Work) ;
		}

		/* ---------------------------------------------------------- */
		/* Numerical & symbolic assembly into current frontal matrix */
		/* ---------------------------------------------------------- */

		UMF_assemble (Numeric, Work) ;

		/* ---------------------------------------------------------- */
		/* scale the pivot column, and save row and column of U and L */
		/* ---------------------------------------------------------- */

		if (!UMF_scale_column (Numeric, Work))
		{
		    return (UMFPACK_ERROR_out_of_memory) ;
		}

		/* ---------------------------------------------------------- */
		/* Numerical update if enough pivots accumulated */
		/* ---------------------------------------------------------- */

		if (Work->fnpiv >= nb)
		{
		    UMF_blas3_update (Work) ;
		}

		Work->pivrow_in_front = FALSE ;
		Work->pivcol_in_front = FALSE ;

		/* ---------------------------------------------------------- */
		/* If front is empty, evaporate it */
		/* ---------------------------------------------------------- */

		if (Work->fnrows == 0 || Work->fncols == 0)
		{
		    /* This does not create an element, just evaporates. */
		    /* It ensures that a front is not 0-by-c or c-by-0. */
		    DEBUG1 (("Evaporate empty front:\n")) ;
		    (void) UMF_create_element (Numeric, Work) ;
		    Work->fnrows = 0 ;
		    Work->fncols = 0 ;
		}
	    }
	}

	/* ------------------------------------------------------------------ */
	/* Wrapup the current frontal matrix.  This is the last in a chain */
	/* in the column elimination tree.  The next frontal matrix */
	/* cannot overlap with the current one, which will be its sibling */
	/* in the column etree. */
	/* ------------------------------------------------------------------ */

	if (!UMF_create_element (Numeric, Work))
	{
	    return (UMFPACK_ERROR_out_of_memory) ;
	}

	/* ------------------------------------------------------------------ */
	/* current front is now empty */
	/* ------------------------------------------------------------------ */

	Work->fnrows = 0 ;
	Work->fncols = 0 ;
    }

    /* ---------------------------------------------------------------------- */
    /* end the last Lchain and Uchain and finalize the LU factors */
    /* ---------------------------------------------------------------------- */

    UMF_kernel_wrapup (Numeric, Symbolic, Work) ;

    /* note that the matrix may be singular */
    return (UMFPACK_OK) ;
}

