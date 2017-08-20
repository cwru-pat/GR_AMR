/* ========================================================================== */
/* === UMF_get_memory ======================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Reallocate the workspace (Numeric->Memory) and shift elements downwards.
    needunits: increase in size so that the free space is at least this many
    Units (to which the tuple lengths is added).

    Return TRUE if successful, FALSE if out of memory.
*/

#include "umf_internal.h"
#include "umf_garbage_collection.h"
#include "umf_tuple_lengths.h"
#include "umf_build_tuples.h"
#include "umf_mem_free_tail_block.h"
#include "umf_realloc.h"

GLOBAL Int UMF_get_memory
(
    NumericType *Numeric,
    WorkType *Work,
    Int needunits
)
{
    Int i, minsize, newsize, newmem, costly, row, col, *Row_tlen, *Col_tlen,
	n_row, n_col, *Row_degree, *Col_degree ;
    Unit *mnew, *p ;
    double nsize, bsize ;

    /* ---------------------------------------------------------------------- */
    /* get and check parameters */
    /* ---------------------------------------------------------------------- */

    ASSERT (Numeric) ;
    ASSERT (Work) ;
    ASSERT (Numeric->Memory) ;

#ifndef NDEBUG
    DEBUG1 (("::::GET MEMORY::::\n")) ;
    UMF_dump_memory (Numeric) ;
#endif

    n_row = Work->n_row ;
    n_col = Work->n_col ;
    Row_degree = Numeric->Rperm ;	/* for NON_PIVOTAL_ROW macro */
    Col_degree = Numeric->Cperm ;	/* for NON_PIVOTAL_COL macro */
    Row_tlen   = Numeric->Uilen ;
    Col_tlen   = Numeric->Lilen ;

    /* ---------------------------------------------------------------------- */
    /* initialize the tuple list lengths */
    /* ---------------------------------------------------------------------- */

    for (row = 0 ; row < n_row ; row++)
    {
	if (NON_PIVOTAL_ROW (row))
	{
	    Row_tlen [row] = 0 ;
	}
    }
    for (col = 0 ; col < n_col ; col++)
    {
	if (NON_PIVOTAL_COL (col))
	{
	    Col_tlen [col] = 0 ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* determine how much memory is needed for the tuples */
    /* ---------------------------------------------------------------------- */

    nsize = (double) needunits + 2 ;
    needunits += UMF_tuple_lengths (Numeric, Work, &nsize) ;
    needunits += 2 ;	/* add 2, so that newmem >= 2 is true if realloc'd */

    /* note: Col_tlen and Row_tlen are updated, but the tuple lists */
    /* themselves are not.  Do not attempt to scan the tuple lists. */
    /* They are now stale, and are about to be destroyed and recreated. */

    /* ---------------------------------------------------------------------- */
    /* determine the desired new size of memory */
    /* ---------------------------------------------------------------------- */

    DEBUG0 (("UMF_get_memory: needunits: "ID"\n", needunits)) ;

    minsize = Numeric->size + needunits ;
    nsize += (double) Numeric->size ;

    bsize = ((double) Int_MAX) / sizeof (Unit) - 1 ;

    if (nsize > bsize)	    /* a double relop, but ignore NaN case */
    {
	return (FALSE) ;	/* the minimum size is larger than Int_MAX */
    }

    newsize = (Int) (UMF_REALLOC_INCREASE * ((double) minsize)) ;
    nsize *= UMF_REALLOC_INCREASE * (1.0 + MAX_EPSILON) ;

    if (newsize < 0 || nsize > bsize)
    {
	newsize = (Int) bsize ;	/* we cannot increase the size beyond bsize */
    }
    else
    {
	newsize = MAX (newsize, minsize) ;
    }

    DEBUG0 ((
    "REALLOC MEMORY: needunits "ID" old size: "ID" new size: "ID" Units \n",
	needunits, Numeric->size, newsize)) ;

    /* Forget where the biggest free block is (we no longer need it) */
    /* since garbage collection will occur shortly. */
    Numeric->ibig = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* reallocate the memory, if possible, and make it bigger */
    /* ---------------------------------------------------------------------- */

    mnew = (Unit *) NULL ;
    while (!mnew)
    {
	mnew = (Unit *) UMF_realloc (Numeric->Memory, newsize, sizeof (Unit)) ;
	if (!mnew)
	{
	    if (newsize == minsize)	/* last realloc attempt failed */
	    {
		/* We failed to get the minimum.  Just stick with the */
		/* current allocation and hope that garbage collection */
		/* can recover enough space. */
		mnew = Numeric->Memory ;	/* no new memory available */
		newsize = Numeric->size ;
	    }
	    else
	    {
		/* otherwise, reduce the request and keep trying */
		newsize = (Int) (UMF_REALLOC_REDUCTION * ((double) newsize)) ;
		newsize = MAX (minsize, newsize) ;
	    }
	}
    }
    ASSERT (mnew) ;

    /* see if realloc had to copy, rather than just extend memory */
    costly = (mnew != Numeric->Memory) ;

    /* ---------------------------------------------------------------------- */
    /* extend the tail portion of memory downwards */
    /* ---------------------------------------------------------------------- */

    Numeric->Memory = mnew ;

    newmem = newsize - Numeric->size ;
    ASSERT (newmem == 0 || newmem >= 2) ;

    if (newmem >= 2)
    {
	/* reallocation succeeded */

	/* point to the old tail marker block of size 1 + header */
	p = Numeric->Memory + Numeric->size - 2 ;

	/* create a new block out of the newly extended memory */
	p->header.size = newmem - 1 ;
	i = Numeric->size - 1 ;
	p += newmem ;

	/* create a new tail marker block */
	p->header.prevsize = newmem - 1 ;
	p->header.size = 1 ;

	Numeric->size = newsize ;

	/* free the new block */
	UMF_mem_free_tail_block (Numeric, i) ;

	Numeric->nrealloc++ ;

	if (costly)
	{
	    Numeric->ncostly++ ;
	}

    }
    DEBUG1 (("Done with realloc memory\n")) ;

    /* ---------------------------------------------------------------------- */
    /* garbage collection on the tail of Numeric->memory (destroys tuples) */
    /* ---------------------------------------------------------------------- */

    UMF_garbage_collection (Numeric, Work) ;

    /* ---------------------------------------------------------------------- */
    /* rebuild the tuples */
    /* ---------------------------------------------------------------------- */

    return (UMF_build_tuples (Numeric, Work)) ;

}

