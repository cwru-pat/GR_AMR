/* ========================================================================== */
/* === umf_internal.h ======================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    This file is for internal use in UMFPACK itself, and should not be included
    in user code.  Use umfpack.h instead.  User-accessible file names and
    routine names all start with the letters "umfpack_".  Non-user-accessible
    file names and routine names all start with "umf_".
*/

/* -------------------------------------------------------------------------- */
/* ANSI standard include files */
/* -------------------------------------------------------------------------- */

/* from stdlib.h:  malloc, free, realloc */
/* and when in debug mode:  rand, RAND_MAX */
#include <stdlib.h>

/* from limits.h:  INT_MAX and LONG_MAX */
#include <limits.h>

/* from float.h:  DBL_EPSILON */
#include <float.h>

/* from stdio.h:  printf, NULL */
/* and when in debug mode:  fopen, fscanf */
#include <stdio.h>

/* from string.h: strcmp */
#include <string.h>

/* from math.h:  sqrt, ceil */
#include <math.h>

/* when debugging, assert.h and the assert macro are used (see umf_dump.h) */

/* -------------------------------------------------------------------------- */
/* MATLAB include files */
/* -------------------------------------------------------------------------- */

#ifdef MATHWORKS
#include "util.h"
#endif

#ifdef MATLAB_MEX_FILE
#include "matrix.h"
#include "mex.h"
#endif

/* -------------------------------------------------------------------------- */
/* Architecture */
/* -------------------------------------------------------------------------- */

#if defined (__sun) || defined (MSOL2) || defined (ARCH_SOL2)
#define UMF_SOL2
#define UMFPACK_ARCHITECTURE "Sun Solaris"

#elif defined (__sgi) || defined (MSGI) || defined (ARCH_SGI)
#define UMF_SGI
#define UMFPACK_ARCHITECTURE "SGI Irix"

#elif defined (__linux) || defined (MGLNX86) || defined (ARCH_GLNX86)
#define UMF_LINUX
#define UMFPACK_ARCHITECTURE "Linux"

#elif defined (_AIX) || defined (MIBM_RS) || defined (ARCH_IBM_RS)
#define UMF_AIX
#define UMFPACK_ARCHITECTURE "IBM AIX"

#elif defined (__alpha) || defined (MALPHA) || defined (ARCH_ALPHA)
#define UMF_ALPHA
#define UMFPACK_ARCHITECTURE "Compaq Alpha"

#elif defined (__WIN32) || defined (_WIN32) || defined (_win32) || defined (__win32) || defined (WIN32)
#define UMF_WINDOWS
#define UMFPACK_ARCHITECTURE "Microsoft Windows"

#elif defined (__hppa) || defined (__hpux) || defined (MHPUX) || defined (ARCH_HPUX)
#define UMF_HP
#define UMFPACK_ARCHITECTURE "HP Unix"

#elif defined (__hp700) || defined (MHP700) || defined (ARCH_HP700)
#define UMF_HP
#define UMFPACK_ARCHITECTURE "HP 700 Unix"

#else
/* If the architecture is unknown, and you call the BLAS, you may need to */
/* define BLAS_BY_VALUE, BLAS_NO_UNDERSCORE, and/or BLAS_CHAR_ARG yourself. */
#define UMFPACK_ARCHITECTURE "unknown"
#endif


/* -------------------------------------------------------------------------- */
/* basic definitions */
/* -------------------------------------------------------------------------- */

#ifdef MAX
#undef MAX
#endif
#ifdef MIN
#undef MIN
#endif

/* for integer MAX/MIN, or for doubles when we don't care how NaN's behave: */
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define ONES_COMPLEMENT(r) (-(r)-1)

/* largest allowable double precision epsilon */
#define MAX_EPSILON 1e-8

/* logical expression of p implies q: */
#define IMPLIES(p,q) (!(p) || (q))

/* Note that the IBM RS 6000 xlc predefines TRUE and FALSE in <types.h>. */
/* The Compaq Alpha also predefines TRUE and FALSE. */
#ifdef TRUE
#undef TRUE
#endif
#ifdef FALSE
#undef FALSE
#endif

#define TRUE (1)
#define FALSE (0)
#define PRIVATE static
#define GLOBAL
#define EMPTY (-1)

/* Note that Linux's gcc 2.96 defines NULL as ((void *) 0), but other */
/* compilers (even gcc 2.95.2 on Solaris) define NULL as 0 or (0). */
#ifdef NULL
#undef NULL
#endif

#define NULL 0

/* -------------------------------------------------------------------------- */
/* Real/complex and int/long definitions, double relops */
/* -------------------------------------------------------------------------- */

#include "umf_version.h"

/* -------------------------------------------------------------------------- */
/* Compile-time configurations */
/* -------------------------------------------------------------------------- */

#include "umf_config.h"

/* -------------------------------------------------------------------------- */
/* umfpack include file */
/* -------------------------------------------------------------------------- */

#include "umfpack.h"

/* -------------------------------------------------------------------------- */
/* for contents of Info.  This must correlate with umfpack.h */
/* -------------------------------------------------------------------------- */

#define ESTIMATE (UMFPACK_NUMERIC_SIZE_ESTIMATE - UMFPACK_NUMERIC_SIZE)
#define ACTUAL 0

/* -------------------------------------------------------------------------- */
/* for clearing the external degree counters */
/* -------------------------------------------------------------------------- */

#define MAX_MARK(n) Int_MAX - (2*(n)+1)

/* -------------------------------------------------------------------------- */
/* dense row/column macro */
/* -------------------------------------------------------------------------- */

/* In order for a row or column to be treated as "dense", it must have more */
/* entries than the value returned by this macro.  n is the dimension of the */
/* matrix, and alpha is the dense row/column control parameter. */

/* Note: this is not defined in alpha is NaN or Inf: */
#define UMFPACK_DENSE_DEGREE_THRESHOLD(alpha,n) \
    ((Int) MAX (16.0, (alpha) * 16.0 * sqrt ((double) (n))))

/* -------------------------------------------------------------------------- */
/* PRINTF */
/* -------------------------------------------------------------------------- */

#define PRINTFk(k,params) { if (prl >= (k)) { PRINTF (params) ; } }
#define PRINTF1(params) PRINTFk (1, params)
#define PRINTF2(params) PRINTFk (2, params)
#define PRINTF3(params) PRINTFk (3, params)
#define PRINTF4(params) PRINTFk (4, params)
#define PRINTF5(params) PRINTFk (5, params)
#define PRINTF6(params) PRINTFk (6, params)

/* -------------------------------------------------------------------------- */
/* Fixed control parameters */
/* -------------------------------------------------------------------------- */

/* maximum number of columns to consider at one time, in a single front */
#define MAX_CANDIDATES 128

/* reduce Numeric->Memory request by this ratio, if allocation fails */
#define UMF_REALLOC_REDUCTION (0.95)

/* increase Numeric->Memory request by this ratio, if we need more */
#define UMF_REALLOC_INCREASE (1.2)

/* -------------------------------------------------------------------------- */
/* Memory space definitions */
/* -------------------------------------------------------------------------- */

/* for memory alignment - assume double has worst case alignment */
typedef double Align ;

/* get number of bytes required to hold n items of a type: */
/* note that this will not overflow, because sizeof (type) is always */
/* greater than or equal to sizeof (Int) >= 2 */
#define BYTES(type,n) (sizeof (type) * (n))

/* ceiling of (b/u).  Assumes b >= 0 and u > 0 */
#define CEILING(b,u) (((b) + (u) - 1) / (u))

/* get number of Units required to hold n items of a type: */
#define UNITS(type,n) (CEILING (BYTES (type, n), sizeof (Unit)))

/* same as DUNITS, but use double instead of int to avoid overflow */
#define DUNITS(type,n) (ceil (BYTES (type, (double) n) / sizeof (Unit)))

union Unit_union
{	/* memory is allocated in multiples of Unit */
    struct
    {
	Int
	    size,	/* size, in Units, of the block, excl. header block */
			/* size >= 0: block is in use */
			/* size < 0: block is free, of |size| Units */
	    prevsize ;	/* size, in Units, of preceding block in S->Memory */
			/* during garbage_collection, prevsize is set to -e-1 */
			/* for element e, or positive (and thus a free block) */
			/* otherwise */
    } header ;		/* block header */
    Align  xxxxxx ;	/* force alignment of blocks (xxxxxx is never used) */
} ;

typedef union Unit_union Unit ;

/* get the size of an allocated block */
#define GET_BLOCK_SIZE(p) (((p)-1)->header.size)

/* -------------------------------------------------------------------------- */
/* Numeric */
/* -------------------------------------------------------------------------- */

/*
    NUMERIC_VALID and SYMBOLIC_VALID:
    The different values of SYBOLIC_VALID and NUMERIC_VALID are chosen as a
    first defense against corrupted *Symbolic or *Numeric pointers passed to an
    UMFPACK routine.  They also ensure that the objects are used only by the
    same version that created them (umfpack_di_*, umfpack_dl_*, umfpack_zi_*,
    or umfpack_zl_*).  The values have also been changed since prior releases of
    the code to ensure that all routines that operate on the objects are of the 
    same release.  The values themselves are purely arbitrary.  The are less
    than the ANSI C required minimums of INT_MAX and LONG_MAX, respectively.
*/

#ifdef DINT
#define NUMERIC_VALID  15674
#define SYMBOLIC_VALID 41234
#endif
#ifdef DLONG
#define NUMERIC_VALID  456789120
#define SYMBOLIC_VALID 432192913
#endif
#ifdef ZINT
#define NUMERIC_VALID  17654
#define SYMBOLIC_VALID 40123
#endif
#ifdef ZLONG
#define NUMERIC_VALID  123987654
#define SYMBOLIC_VALID 119291234
#endif

typedef struct	/* NumericType */
{
    double
	flops,		/* "true" flop count */
	relpt,		/* relative pivot tolerance used */
	relax,		/* relaxed amalgamation parameter */
	relax2,		/* relax2 amalgamation parameter */
	relax3,		/* relax3 amalgamation parameter */
	alloc_init ;	/* initial allocation of Numeric->memory */

    Int valid ;		/* set to NUMERIC_VALID, for validity check */

    /* Memory space for A and LU factors */
    Unit
	*Memory ;	/* working memory for A and LU factors */
    Int
	ihead,		/* pointer to tail of LU factors, in Numeric->Memory */
	itail,		/* pointer to top of elements & tuples,  */
			/* in Numeric->Memory */
	ibig,		/* pointer to largest free block seen in tail */
	size ;		/* size of Memory, in Units */

    Int
	*Rperm,		/* pointer to row perm array, size: n+1 */
			/* after UMF_kernel:  Rperm [new] = old */
			/* during UMF_kernel: Rperm [old] = new */
	*Cperm,		/* pointer to col perm array, size: n+1 */
			/* after UMF_kernel:  Cperm [new] = old */
			/* during UMF_kernel: Cperm [old] = new */

	*Upos,		/* see UMFPACK_get_numeric for a description */
	*Lpos,
	*Lip,
	*Lilen,
	*Uip,
	*Uilen,
	*Upattern ;	/* pattern of last row of U (if singular) */

    Int
	ulen,		/* length of Upattern */
	npiv,		/* number of structural pivots found (sprank approx) */
	nnzpiv ;	/* number of numerical (nonzero) pivots found */

    Entry
	*D ;		/* D [i] is the diagonal entry of U */

    double
	min_udiag,	/* smallest abs value on diagonal of D */
	max_udiag,	/* smallest abs value on diagonal of D */
	rcond ;		/* min (D) / max (D) */

    Int
	n_row, n_col ;	/* A is n_row-by-n_row */

    /* for information only: */
    Int
	tail_usage,	/* amount of memory allocated in tail */
			/* head_usage is Numeric->ihead */
	init_usage,	/* memory usage just after UMF_kernel_init */
	max_usage,	/* peak memory usage (excludes internal and external */
			/* fragmentation in the tail) */
	ngarbage,	/* number of garbage collections performed */
	nrealloc,	/* number of reallocations performed */
	ncostly,	/* number of costly reallocations performed */
	isize,		/* size of integer pattern of L and U */
	nLentries,	/* number of entries in L, excluding diagonal */
	nUentries,	/* number of entries in U, including diagonal */
			/* Some entries may be numerically zero. */
	lnz,		/* number of nonzero entries in L, excl. diagonal */
	unz,		/* number of nonzero entries in U, excl. diagonal */
	maxfrsize ;	/* largest actual front size */

} NumericType ;



/* -------------------------------------------------------------------------- */
/* Element tuples for connecting elements together in a matrix */
/* -------------------------------------------------------------------------- */

typedef struct	/* Tuple */
{
    /* The (e,f) tuples for the element lists */
    Int
	e,		/* element */
	f ;		/* contribution to the row/col appears at this offset */

} Tuple ;

/* Col_degree is aliased with Cperm, and Row_degree with Rperm */
#define NON_PIVOTAL_COL(col) (Col_degree [col] >= 0)
#define NON_PIVOTAL_ROW(row) (Row_degree [row] >= 0)

/* -------------------------------------------------------------------------- */
/* An element */
/* -------------------------------------------------------------------------- */

typedef struct	/* Element */
{
    Int

	cdeg,		/* external column degree + cdeg0 offset */
	rdeg,		/* external row degree    + rdeg0 offset */
	nrowsleft,	/* number of rows remaining */
	ncolsleft,	/* number of columns remaining */
	nrows,		/* number of rows */
	ncols,		/* number of columns */
	next ;		/* for list link of sons, used during assembly only */

    /* followed in memory by:
    Int
	col [0..ncols-1],	column indices of this element
	row [0..nrows-1] ;	row indices of this element
    Entry			(suitably aligned, see macro below)
	C [0...nrows-1, 0...ncols-1] ;
	size of C is nrows*ncols Entry's
    */

} Element ;

/* macros for computing pointers to row/col indices, and contribution block: */

#define GET_ELEMENT_SIZE(nr,nc) \
(UNITS (Element, 1) + UNITS (Int, (nc) + (nr)) + UNITS (Entry, (nc) * (nr)))

#define DGET_ELEMENT_SIZE(nr,nc) \
(DUNITS (Element, 1) + DUNITS (Int, (nc) + (nr)) + DUNITS (Entry, (nc) * (nr)))

#define GET_ELEMENT_COLS(ep,p,Cols) { \
    ASSERT (p) ; \
    ASSERT (p >= Numeric->Memory + Numeric->itail) ; \
    ASSERT (p <= Numeric->Memory + Numeric->size) ; \
    ep = (Element *) p ; \
    p += UNITS (Element, 1) ; \
    Cols = (Int *) p ; \
}

#define GET_ELEMENT_PATTERN(ep,p,Cols,Rows,ncm) { \
    GET_ELEMENT_COLS (ep, p, Cols) ; \
    ncm = ep->ncols ; \
    Rows = Cols + ncm ; \
}

#define GET_ELEMENT(ep,p,Cols,Rows,ncm,nrm,C) { \
    GET_ELEMENT_PATTERN (ep, p, Cols, Rows, ncm) ; \
    nrm = ep->nrows ; \
    p += UNITS (Int, ncm + nrm) ; \
    C = (Entry *) p ; \
}

/* -------------------------------------------------------------------------- */
/* Work data structure */
/* -------------------------------------------------------------------------- */

/*
    This data structure holds items needed only during factorization.
    All of this is freed when UMFPACK_numeric completes.  Note that some of
    it is stored in the tail end of Numeric->S (namely, the Tuples and the
    Elements).
*/

typedef struct	/* WorkType */
{

    /* ---------------------------------------------------------------------- */
    /* information about each row and col of A */
    /* ---------------------------------------------------------------------- */

    /*
	Row_tuples:	pointer to tuple list (alias with Numeric->Uip)
	Row_tlen:	number of tuples (alias with Numeric->Uilen)
	Col_tuples:	pointer to tuple list (alias with Numeric->Lip)
	Col_tlen:	number of tuples (alias with Numeric->Lilen)
	Row_degree:	degree of the row or column (alias Numeric->Rperm)
	Col_degree:	degree of the row or column (alias Numeric->Cperm)

	The Row_degree and Col_degree are MATLAB-style colmmd approximations,
	are equal to the sum of the sizes of the elements (contribution blocks)
	in each row and column.  They are maintained when elements are created
	and assembled.  They are used only during the pivot row and column
	search.  They are not needed to represent the pattern of the remaining
	matrix.
    */

    /* ---------------------------------------------------------------------- */
    /* information about each element */
    /* ---------------------------------------------------------------------- */

    Int	*E ;		/* E [0 .. n_col + n_inner] element "pointers" */
			/* (offsets in Numeric->Memory) */

    /* ---------------------------------------------------------------------- */
    /* generic workspace */
    /* ---------------------------------------------------------------------- */

    Int			/* Sizes:  nn = MAX (n_row, n_col) */
	*Wp,		/* nn+1 */
	*Wm,		/* maxnrows+1 */
	*Wio,		/* maxncols+1 */
	*Woi,		/* maxncols+1 */
	*Woo,		/* MAX (maxnrows,maxncols)+1 */
	*Wcol,		/* pointer to Wm */
	*Wrow ;		/* pointer to Fcols, Wio, or Woi */

    /* ---------------------------------------------------------------------- */
    Int
	*Lpattern,	/* pattern of column of L, for one Lchain */
	*Upattern,	/* pattern of row of U, for one Uchain */
	ulen, llen ;	/* length of Upattern and Lpattern */

    /* ---------------------------------------------------------------------- */

    Int
	n_row, n_col,	/* matrix is n_row-by-n_col */
	nz,		/* nonzeros in the elements for this matrix */
	npiv,		/* number of pivot rows and columns so far */
	Wpflag,
	nel,		/* elements in use are in the range 1..nel */
	nelorig,	/* elements 1..nelorig are original elements */
	prior_element,
	overlap,
	rdeg0, cdeg0,
	rrdeg, ccdeg,
	Candidates [MAX_CANDIDATES],	 /* current candidate pivot columns */
	ncand,		/* number of candidates (some not in Candidates[ ]) */
	nextcand,	/* next candidate to place in Candidate search set */
	pivrow,		/* current pivot row */
	pivcol,		/* current pivot column */
	scan_pivrow, scan_pivcol,
	do_extend,	/* true if the next pivot extends the current front */
	do_update,	/* true if update should be applied */
	do_scan2row,
	do_scan2col,
	pivot_case,
	frontid,	/* id of current frontal matrix */
	nfr,		/* number of frontal matrices */
	maxfrsize,
	maxnrows,	/* largest number of rows in any front */
	maxncols ;	/* largest number of columns in any front */

    /* ---------------------------------------------------------------------- */
    /* For row-merge tree */
    /* ---------------------------------------------------------------------- */

    Int
	*Front_new1strow ;

    /* ---------------------------------------------------------------------- */
    /* current frontal matrix */
    /* ---------------------------------------------------------------------- */

    Entry
	*Fx ;		/* Fx [0..fnrows_max-1, 0..fncols_max-1] */
			/* working array */
    Int
	*Frows,		/* Frows [0.. ]: row indices of F */

	*Fcols,		/* Fcols [0.. ]: column indices of F */

	*Frpos,		/* position of row indices in F, or -1 if not present */
			/* if Frows[i] == row, then Frpos[row] == i */

	*Fcpos,		/* position of col indices in F, or -1 if not present */
			/* if Fcols[j] == col, then */
			/* Fcpos[col] == j*Front->fnrows_max */

	fnrows,		/* number of rows in contribution block in F */
	fncols,		/* number of columns in contribution block in F */
	fnrows_max,	/* column-dimension (max # of rows) of F */
	fncols_max,	/* row-dimension (max # of columns) of F */
	fnpiv,		/* number of pivots in F */
	fcolpos,	/* offset of candidate pivot column in Fcols [...] */
	frowpos,	/* offset of candidate pivot row in Frows [...] */
	fnzeros,	/* number of explicit zero entries in LU block */
	fscan_row,
	fscan_col,
	pivrow_in_front,	/* true if current pivot row in Frows */
	pivcol_in_front ;	/* true if current pivot column in Fcols */

    /* ---------------------------------------------------------------------- */
    /* Current frontal matrix                                                 */
    /* ---------------------------------------------------------------------- */
    /*  Fx points to current frontal matrix (contribution block and LU        */
    /*  factors).  For example, if Front->fnrows = 4, Front->fncols = 6, and  */
    /*  Front->fnpiv = 3, then "x" is a term in the contribution block, "l"   */
    /*  in L1, "u" in U1, "L" in L2, "U" in U2, and "." is unused.  Fx [0] is */
    /*  "X". The first 3 pivot values (diagonal entries in U1) are 1,2, and 3.*/
    /*  For this frontal matrix, Front->fncols_max = 12 (the number of        */
    /*  columns), and Front->fnrows_max = 8 (the number of rows).  The        */
    /*  frontal matrix is (Front->fnrows_max)-by-(Front->fncols_max)          */
    /*									      */
    /*                             |----------- col 1 of L1 and L2, etc.      */
    /*                             V                                          */
    /*       X x x x x x . . . L L L                                          */
    /*       x x x x x x . . . L L L                                          */
    /*       x x x x x x . . . L L L                                          */
    /*       x x x x x x . . . L L L                                          */
    /*       . . . . . . . . . . . .                                          */
    /*       U U U U U U . . . 3 l l         <- row 3 of U1 and U2            */
    /*       U U U U U U . . . u 2 l         <- row 2 of U1 and U2            */
    /*       U U U U U U . . . u u 1         <- row 1 of U1 and U2            */
    /*									      */
    /* ---------------------------------------------------------------------- */

} WorkType ;


/* -------------------------------------------------------------------------- */
/* Symbolic */
/* -------------------------------------------------------------------------- */

/*
    This is is constructed by UMFPACK_symbolic, and is needed by UMFPACK_numeric
    to factor the matrix.
*/

typedef struct	/* SymbolicType */
{

    double
	drow,		/* dense row control parameter used */
	dcol,		/* dense column control parameter used */
	num_mem_usage_est,	/* estimated max Numeric->Memory size */
	num_mem_size_est,	/* estimated final Numeric->Memory size */
	peak_sym_usage ;	/* peak Symbolic and SymbolicWork usage */

    Int valid,		/* set to SYMBOLIC_VALID, for validity check */
	nchains,
	*Chain_start,
	*Chain_maxrows,
	*Chain_maxcols,
	maxfrsize,
	maxnrows,		/* largest number of rows in any front */
	maxncols,		/* largest number of columns in any front */
	*Front_npivcol,		/* Front_npivcol [j] = size of jth supercolumn*/
	*Front_1strow,		/* first row index in front j */
	*Front_leftmostdesc,	/* leftmost desc of front j */
	*Front_parent,		/* super-column elimination tree */
	*Cperm_init,		/* initial column ordering */
	*Rperm_init,		/* initial row ordering */
	nfr,
	n_row, n_col,		/* matrix A is n_row-by-n_col */
	nz,			/* nz of original matrix */
	nb,			/* block size for BLAS 3 */
	num_mem_init_usage ;	/* min Numeric->Memory for UMF_kernel_init */

} SymbolicType ;


/* -------------------------------------------------------------------------- */
/* for debugging only: */
/* -------------------------------------------------------------------------- */

#include "umf_dump.h"

