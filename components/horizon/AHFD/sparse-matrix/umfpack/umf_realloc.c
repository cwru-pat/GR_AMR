/* ========================================================================== */
/* === UMF_realloc ========================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Realloc a block previously allocated by UMF_malloc.
    Return NULL on failure (in which case the block is still allocated, and will
    be kept at is present size).  This routine is only used for Numeric->Memory.
*/

#include "umf_internal.h"



GLOBAL void *UMF_realloc
(
    void *p,
    Int n_objects,
    size_t size_of_object
)
{
    size_t size ;

#ifdef UMF_TCOV_TEST
    /* For exhaustive statement coverage testing only! */
    /* Pretend to fail, to test out-of-memory conditions. */
    umf_realloc_fail-- ;
    if (umf_realloc_fail <= umf_realloc_hi && umf_realloc_fail >= umf_realloc_lo) { return ((void *) NULL) ; }
#endif

    /* make sure that we allocate something */
    n_objects = MAX (1, n_objects) ;

    size = (size_t) n_objects ;
    ASSERT (size_of_object > 1) ;
    if (size > Int_MAX / size_of_object)
    {
	/* object is too big for integer pointer arithmetic */
	return ((void *) NULL) ;
    }
    size *= size_of_object ;

    DEBUG0 (("UMF_realloc n_objects "ID"  size_of_object "ID"\n",
	n_objects, (Int) size_of_object)) ;

    /* see umf_config.h for the memory allocator selection */
    return (REALLOCATE (p, size)) ;
}

