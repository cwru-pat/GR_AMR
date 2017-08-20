/* ========================================================================== */
/* === UMF_malloc =========================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Allocate a block of n objects, each of a given size.  This routine does not
    handle the case when the size is 1 (allocating char's) because of potential
    integer overflow.  UMFPACK never does that.
    Also maintains the UMFPACK malloc count.
*/

#include "umf_internal.h"

#if defined (UMF_MALLOC_COUNT) || !defined (NDEBUG)

/*
    UMF_malloc_count is a count of the objects malloc'd by UMFPACK.
    It is increased by 7 by UMFPACK_*symbolic, and by 15 or 16 by
    UMFPACK_*numeric.  It is reduced by the same amount by the corresponding
    UMFPACK_free_* routines.  If you suspect a memory leak in your program
    (caused by not properly destroying the Symbolic and Numeric objects)
    then compile with -DUMF_MALLOC_COUNT and check value of UMF_malloc_count.
    By default, UMF_MALLOC_COUNT is not defined, and thus UMFPACK has
    no global variables.
*/

GLOBAL Int UMF_malloc_count = 0 ;

#endif


GLOBAL void *UMF_malloc
(
    Int n_objects,
    size_t size_of_object
)
{
    size_t size ;
    void *p ;

#ifdef UMF_TCOV_TEST
    /* For exhaustive statement coverage testing only! */
    /* Pretend to fail, to test out-of-memory conditions. */
    umf_fail-- ;
    if (umf_fail <= umf_fail_hi && umf_fail >= umf_fail_lo) { return ((void *) NULL) ; }
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

    /* see umf_config.h for the memory allocator selection */
    p = ALLOCATE (size) ;

#if defined (UMF_MALLOC_COUNT) || !defined (NDEBUG)
    if (p)
    {
	/* One more object has been malloc'ed.  Keep track of the count. */
	/* (purely for sanity checks). */
	UMF_malloc_count++ ;
    }
#endif

    return (p) ;
}

