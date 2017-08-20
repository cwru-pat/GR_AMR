/* ========================================================================== */
/* === UMFPACK_free_numeric ================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*  User-callable.  Free the entire Numeric object.
    See UMFPACK_free_numeric.h for details.
*/

#include "umf_internal.h"
#include "umf_free.h"

GLOBAL void UMFPACK_free_numeric
(
    void **NumericHandle
)
{

    NumericType *Numeric ;
    if (!NumericHandle)
    {
	return ;
    }
    Numeric = *((NumericType **) NumericHandle) ;
    if (!Numeric)
    {
	return ;
    }

    (void) UMF_free ((void *) Numeric->D) ;
    (void) UMF_free ((void *) Numeric->Rperm) ;
    (void) UMF_free ((void *) Numeric->Cperm) ;
    (void) UMF_free ((void *) Numeric->Lpos) ;
    (void) UMF_free ((void *) Numeric->Lilen) ;
    (void) UMF_free ((void *) Numeric->Lip) ;
    (void) UMF_free ((void *) Numeric->Upos) ;
    (void) UMF_free ((void *) Numeric->Uilen) ;
    (void) UMF_free ((void *) Numeric->Uip) ;
    (void) UMF_free ((void *) Numeric->Upattern) ;
    (void) UMF_free ((void *) Numeric->Memory) ;
    (void) UMF_free ((void *) Numeric) ;

    *NumericHandle = (void *) NULL ;
}

