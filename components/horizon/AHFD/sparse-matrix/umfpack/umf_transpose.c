/* ========================================================================== */
/* === UMF_transpose ======================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*  Not user-callable.  Computes a permuted transpose, R = (A (P,Q(1:nq)))' in
	MATLAB notation, where R is in column-form.  A is n_row-by-n_col, the
	row-form matrix R is n_row-by-nq, where nq <= n_col.  A may be singular.
	The complex version can do transpose (') or array transpose (.').

	Uses Gustavson's method (Two Fast Algorithms for Sparse Matrices:
	Multiplication and Permuted Transposition, ACM Trans. on Math. Softw.,
	vol 4, no 3, pp. 250-269).
*/

#include "umf_internal.h"
#include "umf_is_permutation.h"

/* ========================================================================== */

GLOBAL Int UMF_transpose
(
    Int n_row,			/* A is n_row-by-n_col */
    Int n_col,
    const Int Ap [ ],		/* size n_col+1 */
    const Int Ai [ ],		/* size nz = Ap [n_col] */
    const double Ax [ ],	/* size nz if present */

    const Int P [ ],	/* P [k] = i means original row i is kth row in A(P,Q)*/
			/* P is identity if not present */
			/* size n_row, if present */

    const Int Q [ ],	/* Q [k] = j means original col j is kth col in A(P,Q)*/
			/* Q is identity if not present */
			/* size nq, if present */
    Int nq,		/* size of Q, ignored if Q is (Int *) NULL */

			/* output matrix: Rp, Ri, Rx, and Rz: */
    Int Rp [ ],		/* size n_row+1 */
    Int Ri [ ],		/* size nz */
    double Rx [ ],	/* size nz, if present */

    Int W [ ],		/* size max (n_row,n_col) workspace */

    Int check		/* if true, then check inputs */
#ifdef COMPLEX
    , const double Az [ ]	/* size nz */
    , double Rz [ ]		/* size nz */
    , Int do_conjugate		/* if true, then do conjugate transpose */
				/* otherwise, do array transpose */
#endif
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int i, j, k, p, bp, nz, newj, ilast, do_values ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    ASSERT (n_col >= 0) ;
    nz = Ap ? Ap [n_col] : 0 ;
    DEBUG2 (("UMF_transpose:  nz "ID"\n", nz)) ;
    DEBUG2 (("n_row "ID"  & "ID"\n", n_row, &n_row)) ;
    DEBUG2 (("n_col "ID"  & "ID"\n", n_col, &n_col)) ;
    DEBUG2 (("Ap   & "ID" to & "ID"\n", Ap, Ap + (n_col+1) - 1)) ;
    DEBUG2 (("Ai   & "ID" to & "ID"\n", Ai, Ai + (nz) - 1)) ;
    if (Ax) DEBUG2 (("Ax   & "ID" to & "ID"\n", Ax, Ax + (nz) - 1)) ;
    if (P)  DEBUG2 (("P    & "ID" to & "ID"\n", P , P  + (n_row) - 1)) ;
    if (Q)  DEBUG2 (("Q    & "ID" to & "ID"\n", Q , Q  + (n_col) - 1)) ;
    DEBUG2 (("Rp   & "ID" to & "ID"\n", Rp, Rp + (n_row+1) - 1)) ;
    DEBUG2 (("Ri   & "ID" to & "ID"\n", Ri, Ri + (nz) - 1)) ;
    if (Rx) DEBUG2 (("Rx   & "ID" to & "ID"\n", Rx, Rx + (nz) - 1)) ;
    DEBUG2 (("W    & "ID" to & "ID"\n", W, W + MAX(n_row,n_col) - 1)) ;
#ifdef COMPLEX
    if (Az) DEBUG2 (("Az   & "ID" to & "ID"\n", Az, Az + (nz) - 1)) ;
    if (Rz) DEBUG2 (("Rz   & "ID" to & "ID"\n", Rz, Rz + (nz) - 1)) ;
    DEBUG2 (("do_conjugate "ID"  & "ID"\n", do_conjugate, &do_conjugate)) ;
#endif
#endif

    if (check)
    {
	/* UMFPACK_symbolic skips this check */
	/* UMFPACK_transpose always does this check */

	if (!Ai || !Ap || !Ri || !Rp || !W)
	{
	    return (UMFPACK_ERROR_argument_missing) ;
	}

	if (n_row <= 0 || n_col <= 0)		/* n_row,n_col must be > 0 */
	{
	    return (UMFPACK_ERROR_n_nonpositive) ;
	}

	nz = Ap [n_col] ;
	if (nz < 0)		/* nz must be >= 0 */
	{
	    return (UMFPACK_ERROR_nz_negative) ;
	}

	if (!UMF_is_permutation (P, W, n_row, n_row) ||
	    !UMF_is_permutation (Q, W, nq, nq))
	{
	    return (UMFPACK_ERROR_invalid_permutation) ;
	}

	if (Ap [0] != 0)
	{
	    return (UMFPACK_ERROR_Ap0_nonzero) ;
	}

	for (j = 0 ; j < n_col ; j++)
	{
	    if (Ap [j] > Ap [j+1])
	    {
		return (UMFPACK_ERROR_col_length_negative) ;
	    }
	}

	for (j = 0 ; j < n_col ; j++)
	{
	    ilast = -1 ;
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		i = Ai [p] ;
		if (i < 0 || i >= n_row)
		{
		    return (UMFPACK_ERROR_row_index_out_of_bounds) ;
		}
		if (i <= ilast)
		{
		    return (UMFPACK_ERROR_jumbled_matrix) ;
		}
		ilast = i ;
	    }
	}
    }

#ifndef NDEBUG
    DEBUG2 (("UMF_transpose, input matrix:\n")) ;
    UMF_dump_col_matrix (Ax,
#ifdef COMPLEX
	Az,
#endif
	Ai, Ap, n_row, n_col, nz) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* count the entries in each row of A */
    /* ---------------------------------------------------------------------- */

    /* use W as workspace for RowCount */

    for (i = 0 ; i < n_row ; i++)
    {
	W [i] = 0 ;
	Rp [i] = 0 ;
    }

    if (Q)
    {
        for (newj = 0 ; newj < nq ; newj++)
        {
	    j = Q [newj] ;
	    ASSERT (j >= 0 && j < n_col) ;
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		i = Ai [p] ;
		ASSERT (i >= 0 && i < n_row) ;
	        W [i]++ ;
	    }
        }
    }
    else
    {
        for (j = 0 ; j < n_col ; j++)
        {
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		i = Ai [p] ;
		ASSERT (i >= 0 && i < n_row) ;
	        W [i]++ ;
	    }
        }
    }

    /* ---------------------------------------------------------------------- */
    /* compute the row pointers for R = A (P,Q) */
    /* ---------------------------------------------------------------------- */

    if (P)
    {
	Rp [0] = 0 ;
	for (k = 0 ; k < n_row ; k++)
	{
	    i = P [k] ;
	    ASSERT (i >= 0 && i < n_row) ;
	    Rp [k+1] = Rp [k] + W [i] ;
	}
	for (k = 0 ; k < n_row ; k++)
	{
	    i = P [k] ;
	    ASSERT (i >= 0 && i < n_row) ;
	    W [i] = Rp [k] ;
	}
    }
    else
    {
	Rp [0] = 0 ;
	for (i = 0 ; i < n_row ; i++)
	{
	    Rp [i+1] = Rp [i] + W [i] ;
	}
	for (i = 0 ; i < n_row ; i++)
	{
	    W [i] = Rp [i] ;
	}
    }
    ASSERT (Rp [n_row] <= Ap [n_col]) ;

    /* at this point, W holds the permuted row pointers */

    /* ---------------------------------------------------------------------- */
    /* construct the row form of B */
    /* ---------------------------------------------------------------------- */

    do_values = Ax && Rx ;
#ifdef COMPLEX
    do_values = do_values && Az && Rz ;
#endif

#ifdef COMPLEX
    if (do_conjugate && do_values)
    {
	if (Q)
	{

		/* R = A (P,Q)' */
		for (newj = 0 ; newj < nq ; newj++)
		{
		    j = Q [newj] ;
		    ASSERT (j >= 0 && j < n_col) ;
		    for (p = Ap [j] ; p < Ap [j+1] ; p++)
		    {
			bp = W [Ai [p]]++ ;
			Ri [bp] = newj ;
			Rx [bp] = Ax [p] ;
			Rz [bp] = -Az [p] ;
		    }
		}

	}
	else
	{

		/* R = A (P,:)' */
		for (j = 0 ; j < n_col ; j++)
		{
		    for (p = Ap [j] ; p < Ap [j+1] ; p++)
		    {
			bp = W [Ai [p]]++ ;
			Ri [bp] = j ;
			Rx [bp] = Ax [p] ;
			Rz [bp] = -Az [p] ;
		    }
		}

	}
    }
    else
#endif
    {
	if (Q)
	{
	    if (do_values)
	    {

		/* R = A (P,Q).' */
		for (newj = 0 ; newj < nq ; newj++)
		{
		    j = Q [newj] ;
		    ASSERT (j >= 0 && j < n_col) ;
		    for (p = Ap [j] ; p < Ap [j+1] ; p++)
		    {
			bp = W [Ai [p]]++ ;
			Ri [bp] = newj ;
			Rx [bp] = Ax [p] ;
#ifdef COMPLEX
			Rz [bp] = Az [p] ;
#endif
		    }
		}

	    }
	    else
	    {

		/* R = pattern of A (P,Q).' */
		for (newj = 0 ; newj < nq ; newj++)
		{
		    j = Q [newj] ;
		    ASSERT (j >= 0 && j < n_col) ;
		    for (p = Ap [j] ; p < Ap [j+1] ; p++)
		    {
			Ri [W [Ai [p]]++] = newj ;
		    }
		}

	    }
	}
	else
	{
	    if (do_values)
	    {

		/* R = A (P,:).' */
		for (j = 0 ; j < n_col ; j++)
		{
		    for (p = Ap [j] ; p < Ap [j+1] ; p++)
		    {
			bp = W [Ai [p]]++ ;
			Ri [bp] = j ;
			Rx [bp] = Ax [p] ;
#ifdef COMPLEX
			Rz [bp] = Az [p] ;
#endif
		    }
		}

	    }
	    else
	    {

		/* R = pattern of A (P,:).' */
		for (j = 0 ; j < n_col ; j++)
		{
		    for (p = Ap [j] ; p < Ap [j+1] ; p++)
		    {
			Ri [W [Ai [p]]++] = j ;
		    }
		}

	    }
	}

    }

#ifndef NDEBUG
    for (k = 0 ; k < n_row ; k++)
    {
	if (P)
	{
	    i = P [k] ;
	}
	else
	{
	    i = k ;
	}
	DEBUG3 ((ID":  W[i] "ID" Rp[k+1] "ID"\n", i, W [i], Rp [k+1])) ;
	ASSERT (W [i] == Rp [k+1]) ;
    }
    DEBUG2 (("UMF_transpose, output matrix:\n")) ;
    UMF_dump_col_matrix (Rx,
#ifdef COMPLEX
	Rz,
#endif
	Ri, Rp, n_col, n_row, Rp [n_row]) ;
#endif

    return (UMFPACK_OK) ;
}

