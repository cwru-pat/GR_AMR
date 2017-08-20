/* ========================================================================== */
/* === UMF_free ============================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Free a block previously allocated by UMF_malloc and return NULL.
    Usage is p = UMF_free (p), to ensure that we don't free it twice.
    Also maintains the UMFPACK malloc count.
*/

#include "umf_internal.h"

#if defined (UMF_MALLOC_COUNT) || !defined (NDEBUG)
#include "umf_malloc.h"
#endif

GLOBAL void *UMF_free
(
    void *p
)
{
    if (p)
    {

	/* see umf_config.h for the memory allocator selection */
	FREE (p) ;

#if defined (UMF_MALLOC_COUNT) || !defined (NDEBUG)
	/* One more object has been free'd.  Keep track of the count. */
	/* (purely for sanity checks). */
	UMF_malloc_count-- ;
#endif

    }

    return ((void *) NULL) ;
}

