/* ========================================================================== */
/* === UMF_order_front_tree ================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Post-ordering of supernodal column elimination tree.
*/

#include "umf_internal.h"

GLOBAL Int UMF_order_front_tree
(
    Int root,
    Int k,
    Int Front_child [ ],		/* input argument, destroyed */
    const Int Front_sibling [ ],
    Int Front_order [ ],
    Int Stack [ ]
)
{
    Int f, head, h, i ;

/* recursive version (Stack [ ] is not used):
    i = root ;
    for (f = Front_child [i] ; f != EMPTY ; f = Front_sibling [f])
    {
	k = UMF_order_front_tree (f, k, Front_child, Front_sibling,
	    Front_order) ;
    }
    Front_order [i] = k++ ;
    return (k) ;
*/

    /* push root on the stack */
    head = 0 ;
    Stack [0] = root ;

    while (head >= 0)
    {
	/* get head of stack */
	i = Stack [head] ;
	DEBUG1 (("head of stack "ID" \n", i)) ;
	ASSERT (i >= 0 && i < UMF_nbug && head <= UMF_fbug) ;

	if (Front_child [i] != EMPTY)
	{
	    /* the children of i are not yet ordered */
	    /* push each child onto the stack in reverse order */
	    /* so that small ones at the head of the list get popped first */
	    /* and the biggest one at the end of the list gets popped last */
	    for (f = Front_child [i] ; f != EMPTY ; f = Front_sibling [f])
	    {
		head++ ;
	    }
	    h = head ;
	    ASSERT (head <= UMF_fbug) ;
	    for (f = Front_child [i] ; f != EMPTY ; f = Front_sibling [f])
	    {
		Stack [h--] = f ;
		DEBUG1 (("push "ID" on stack\n", f)) ;
		ASSERT (f >= 0 && f < UMF_nbug) ;
	    }
	    ASSERT (Stack [h] == i) ;

	    /* delete child list so that i gets ordered next time we see it */
	    Front_child [i] = EMPTY ;
	}
	else
	{
	    /* the children of i (if there were any) are already ordered */
	    /* remove i from the stack and order it.  Front i is kth front */
	    head-- ;
	    DEBUG1 (("pop "ID" order "ID"\n", i, k)) ;
	    Front_order [i] = k++ ;
	    ASSERT (k <= UMF_fbug) ;
	}

#ifndef NDEBUG
	DEBUG1 (("\nStack:")) ;
	for (h = head ; h >= 0 ; h--)
	{
	    Int j = Stack [h] ;
	    DEBUG1 ((" "ID, j)) ;
	    ASSERT (j >= 0 && j < UMF_nbug) ;
	}
	DEBUG1 (("\n\n")) ;
	ASSERT (head < UMF_fbug) ;
#endif

    }
    return (k) ;
}

