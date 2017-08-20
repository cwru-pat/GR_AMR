/* ========================================================================== */
/* === umf_version ========================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
   Define routine names and integer type, depending on version being compiled.

   DINT:	double precision, int's as integers
   DLONG:	double precision, long's as integers

   ZLONG:	complex double precision, long's as integers
   ZINT:	complex double precision, int's as integers
*/

/* Set DINT as the default, if nothing is defined */
#if !defined (DLONG) && !defined (DINT) && !defined (ZLONG) && !defined (ZINT)
#define DINT
#endif

/* Determine if this is a real or complex version */
#if defined (ZLONG) || defined (ZINT)
#define COMPLEX
#endif

/* -------------------------------------------------------------------------- */
/* integer type */
/* -------------------------------------------------------------------------- */

#if defined (DLONG) || defined (ZLONG)

#define LONG_INTEGER
#define Int long
#define ID "%ld"
#define Int_MAX LONG_MAX

#else

#define Int int
#define ID "%d"
#define Int_MAX INT_MAX

#endif

/* -------------------------------------------------------------------------- */
/* Numerical relop macros for correctly handling the NaN case */
/* -------------------------------------------------------------------------- */

/*
SCALAR_IS_NAN(x):
    True if x is NaN.  False otherwise.  The commonly-existing isnan(x)
    function could be used, but it's not in Kernighan & Ritchie 2nd edition
    (ANSI C).  It may appear in <math.h>, but I'm not certain about
    portability.  The expression x != x is true if and only if x is NaN,
    according to the IEEE 754 floating-point standard.

SCALAR_IS_ZERO(x):
    True if x is zero.  False if x is nonzero, NaN, or +/- Inf.
    This is (x == 0) if the compiler is IEEE 754 compliant.

SCALAR_IS_NONZERO(x):
    True if x is nonzero, NaN, or +/- Inf.  False if x zero.
    This is (x != 0) if the compiler is IEEE 754 compliant.

SCALAR_IS_LTZERO(x):
    True if x is < zero or -Inf.  False if x is >= 0, NaN, or +Inf.
    This is (x < 0) if the compiler is IEEE 754 compliant.
*/

#if defined (MATHWORKS)

/* The MathWorks has their own macros in util.h that handle NaN's properly. */
#define SCALAR_IS_NAN(x)	(utIsNaN (x))
#define SCALAR_IS_ZERO(x)	(utEQZero (x))
#define SCALAR_IS_NONZERO(x)	(utNEZero (x))
#define SCALAR_IS_LTZERO(x)	(utLTZero (x))

#elif defined (UMF_WINDOWS)

/* Yes, this is exceedingly ugly.  Blame Microsoft, which hopelessly */
/* violates the IEEE 754 floating-point standard in a bizarre way. */
/* If you're using an IEEE 754-compliant compiler, then x != x is true */
/* iff x is NaN.  For Microsoft, (x < x) is true iff x is NaN. */
/* So either way, this macro safely detects a NaN. */
#define SCALAR_IS_NAN(x)	(((x) != (x)) || (((x) < (x))))
#define SCALAR_IS_ZERO(x)	(((x) == 0.) && !SCALAR_IS_NAN(x))
#define SCALAR_IS_NONZERO(x)	(((x) != 0.) || SCALAR_IS_NAN(x))
#define SCALAR_IS_LTZERO(x)	(((x) < 0.) && !SCALAR_IS_NAN(x))

#else

/* These all work properly, according to the IEEE 754 standard ... except on */
/* a PC with windows.  Works fine in Linux on the same PC... */
#define SCALAR_IS_NAN(x)	((x) != (x))
#define SCALAR_IS_ZERO(x)	((x) == 0.)
#define SCALAR_IS_NONZERO(x)	((x) != 0.)
#define SCALAR_IS_LTZERO(x)	((x) < 0.)

#endif

/* scalar absolute value macro. If x is NaN, the result is NaN: */
#define SCALAR_ABS(x) ((SCALAR_IS_LTZERO (x)) ? -(x) : (x))

/* true if an integer (stored in double x) would overflow (or if x is NaN) */
#define INT_OVERFLOW(x) ((!((x) * (1.0 + MAX_EPSILON) <= Int_MAX)) || SCALAR_IS_NAN (x))


/* -------------------------------------------------------------------------- */
/* Real floating-point arithmetic */
/* -------------------------------------------------------------------------- */

#ifndef COMPLEX

#define Entry double

#define REAL_COMPONENT(c) 	(c)
#define IMAG_COMPONENT(c)	(0.)
#define ASSIGN(c,s1,s2)		{ (c) = (s1) ; }
#define CLEAR(c)		{ (c) = 0. ; }
#define IS_ZERO(a)		SCALAR_IS_ZERO (a)
#define IS_NONZERO(a)		SCALAR_IS_NONZERO (a)
#define ASSEMBLE(c,a)		{ (c) += (a) ; }
#define DECREMENT(c,a)		{ (c) -= (a) ; }
#define MULT(c,a,b)		{ (c) = (a) * (b) ; }
#define MULT_CONJ(c,a,b)	{ (c) = (a) * (b) ; }
#define MULT_SUB(c,a,b)		{ (c) -= (a) * (b) ; }
#define MULT_SUB_CONJ(c,a,b)	{ (c) -= (a) * (b) ; }
#define DIV(c,a,b)		{ (c) = (a) / (b) ; }
#define DIV_CONJ(c,a,b)		{ (c) = (a) / (b) ; }
#define APPROX_ABS(s,a)		{ (s) = SCALAR_ABS (a) ; }
#define ABS(s,a)		{ (s) = SCALAR_ABS (a) ; }

/* print an entry (avoid printing "-0" for negative zero).  */
#define PRINT_ENTRY(a) \
{ \
    if (SCALAR_IS_NONZERO (a)) \
    { \
	PRINTF ((" (%g)", (a))) ; \
    } \
    else \
    { \
	PRINTF ((" (0)")) ; \
    } \
}

/* for flop counts */
#define MULTSUB_FLOPS	2.	/* c -= a*b */
#define DIV_FLOPS	1.	/* c = a/b */
#define ABS_FLOPS	0.	/* c = abs (a) */
#define ASSEMBLE_FLOPS	1.	/* c += a */
#define DECREMENT_FLOPS	1.	/* c -= a */
#define MULT_FLOPS	1.	/* c = a*b */

#else

/* -------------------------------------------------------------------------- */
/* Complex floating-point arithmetic */
/* -------------------------------------------------------------------------- */

/*
    Note:  An alternative to this DoubleComplex type would be to use a
    struct { double r ; double i ; }.  The problem with that method
    (used by the Sun Performance Library, for example) is that ANSI C provides
    no guarantee about the layout of a struct.  It is possible that the sizeof
    the struct above would be greater than 2 * sizeof (double).  This would
    mean that the complex BLAS could not be used.  The method used here avoids
    that possibility.  ANSI C *does* guarantee that an array of structs has
    the same size as n times the size of one struct.

    The ANSI C99 version of the C language includes a "double _Complex" type.
    It should be possible in that case to do the following:

    #define Entry double _Complex

    and remove the DoubleComplex struct.  The macros, below, could then be
    replaced with instrinsic operators.  Note that the #define Real and
    #define Imag should also be removed (they only appear in this file).

    For the MULT, MULT_SUB, MULT_SUB_CONJ, DIV, and DIV_CONJ macros,
    the output argument c cannot be the same as any input argument.

*/

typedef struct
{
    double component [2] ;	/* real and imaginary parts */

} DoubleComplex ;

#define Entry DoubleComplex
#define Real component [0]
#define Imag component [1]

/* for flop counts */
#define MULTSUB_FLOPS	8.	/* c -= a*b */
#define DIV_FLOPS	9.	/* c = a/b */
#define ABS_FLOPS	6.	/* c = abs (a), count sqrt as one flop */
#define ASSEMBLE_FLOPS	2.	/* c += a */
#define DECREMENT_FLOPS	2.	/* c -= a */
#define MULT_FLOPS	6.	/* c = a*b */

/* -------------------------------------------------------------------------- */

/* real part of c */
#define REAL_COMPONENT(c) ((c).Real)

/* -------------------------------------------------------------------------- */

/* imag part of c */
#define IMAG_COMPONENT(c) ((c).Imag)

/* -------------------------------------------------------------------------- */

/* c = (s1) + (s2)i */
#define ASSIGN(c,s1,s2) \
{ \
    (c).Real = (s1) ; \
    (c).Imag = (s2) ; \
}

/* -------------------------------------------------------------------------- */

/* c = 0 */
#define CLEAR(c) \
{ \
    (c).Real = 0. ; \
    (c).Imag = 0. ; \
}

/* -------------------------------------------------------------------------- */

/* True if a == 0 */
#define IS_ZERO(a) \
    (SCALAR_IS_ZERO ((a).Real) && SCALAR_IS_ZERO ((a).Imag))

/* -------------------------------------------------------------------------- */

/* True if a != 0 */
#define IS_NONZERO(a) \
    (SCALAR_IS_NONZERO ((a).Real) || SCALAR_IS_NONZERO ((a).Imag))

/* -------------------------------------------------------------------------- */

/* c += a */
#define ASSEMBLE(c,a) \
{ \
    (c).Real += (a).Real ; \
    (c).Imag += (a).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c -= a */
#define DECREMENT(c,a) \
{ \
    (c).Real -= (a).Real ; \
    (c).Imag -= (a).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c = a*b, assert because c cannot be the same as a or b */
#define MULT(c,a,b) \
{ \
    ASSERT (&(c) != &(a) && &(c) != &(b)) ; \
    (c).Real = (a).Real * (b).Real - (a).Imag * (b).Imag ; \
    (c).Imag = (a).Imag * (b).Real + (a).Real * (b).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c = a*conjugate(b), assert because c cannot be the same as a or b */
#define MULT_CONJ(c,a,b) \
{ \
    ASSERT (&(c) != &(a) && &(c) != &(b)) ; \
    (c).Real = (a).Real * (b).Real + (a).Imag * (b).Imag ; \
    (c).Imag = (a).Imag * (b).Real - (a).Real * (b).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c -= a*b, assert because c cannot be the same as a or b */
#define MULT_SUB(c,a,b) \
{ \
    ASSERT (&(c) != &(a) && &(c) != &(b)) ; \
    (c).Real -= (a).Real * (b).Real - (a).Imag * (b).Imag ; \
    (c).Imag -= (a).Imag * (b).Real + (a).Real * (b).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c -= a*conjugate(b), assert because c cannot be the same as a or b */
#define MULT_SUB_CONJ(c,a,b) \
{ \
    ASSERT (&(c) != &(a) && &(c) != &(b)) ; \
    (c).Real -= (a).Real * (b).Real + (a).Imag * (b).Imag ; \
    (c).Imag -= (a).Imag * (b).Real - (a).Real * (b).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c = a/b, be careful to avoid underflow and overflow */
#ifdef MATHWORKS
#define DIV(c,a,b) \
{ \
    (void) utDivideComplex ((a).Real, (a).Imag, (b).Real, (b).Imag, \
	&((c).Real), &((c).Imag)) ; \
}
#else
/* This uses ACM Algo 116, by R. L. Smith, 1962. */
/* c cannot be the same variable as a or b. */
/* Ignore NaN case for double relop br>=bi. */
#define DIV(c,a,b) \
{ \
    double r, den, br, bi ; \
    ASSERT (&(c) != &(a) && &(c) != &(b)) ; \
    br = SCALAR_ABS ((b).Real) ; \
    bi = SCALAR_ABS ((b).Imag) ; \
    if (br >= bi) \
    { \
        r = (b).Imag / (b).Real ; \
        den = (b).Real + r * (b).Imag ; \
        (c).Real = ((a).Real + (a).Imag * r) / den ; \
        (c).Imag = ((a).Imag - (a).Real * r) / den ; \
    } \
    else \
    { \
        r = (b).Real / (b).Imag ; \
        den = r * (b).Real + (b).Imag ; \
        (c).Real = ((a).Real * r + (a).Imag) / den ; \
        (c).Imag = ((a).Imag * r - (a).Real) / den ; \
    } \
}
#endif

/* -------------------------------------------------------------------------- */

/* c = a/conjugate(b), be careful to avoid underflow and overflow */
#ifdef MATHWORKS
#define DIV_CONJ(c,a,b) \
{ \
    (void) utDivideComplex ((a).Real, (a).Imag, (b).Real, (-(b).Imag), \
	&((c).Real), &((c).Imag)) ; \
}
#else
/* This uses ACM Algo 116, by R. L. Smith, 1962. */
/* c cannot be the same variable as a or b. */
/* Ignore NaN case for double relop br>=bi. */
#define DIV_CONJ(c,a,b) \
{ \
    double r, den, br, bi ; \
    ASSERT (&(c) != &(a) && &(c) != &(b)) ; \
    br = SCALAR_ABS ((b).Real) ; \
    bi = SCALAR_ABS ((b).Imag) ; \
    if (br >= bi) \
    { \
        r = (-(b).Imag) / (b).Real ; \
        den = (b).Real - r * (b).Imag ; \
        (c).Real = ((a).Real + (a).Imag * r) / den ; \
        (c).Imag = ((a).Imag - (a).Real * r) / den ; \
    } \
    else \
    { \
        r = (b).Real / (-(b).Imag) ; \
        den =  r * (b).Real - (b).Imag; \
        (c).Real = ((a).Real * r + (a).Imag) / den ; \
        (c).Imag = ((a).Imag * r - (a).Real) / den ; \
    } \
}
#endif

/* -------------------------------------------------------------------------- */

/* approximate absolute value, s = |r|+|i| */
#define APPROX_ABS(s,a) \
{ \
    (s) = SCALAR_ABS ((a).Real) + SCALAR_ABS ((a).Imag) ; \
}

/* -------------------------------------------------------------------------- */

/* exact absolute value, s = sqrt (a.real^2 + amag^2) */
#ifdef MATHWORKS
#define ABS(s,a) \
{ \
    (s) = utFdlibm_hypot ((a).Real, (a).Imag) ; \
}
#else
/* Ignore NaN case for the double relops ar>=ai and ar+ai==ar. */
#define ABS(s,a) \
{ \
    double r, ar, ai ; \
    ar = SCALAR_ABS ((a).Real) ; \
    ai = SCALAR_ABS ((a).Imag) ; \
    if (ar >= ai) \
    { \
	if (ar + ai == ar) \
	{ \
	    (s) = ar ; \
	} \
	else \
	{ \
	    r = ai / ar ; \
	    (s) = ar * sqrt (1.0 + r*r) ; \
	} \
    } \
    else \
    { \
	if (ai + ar == ai) \
	{ \
	    (s) = ai ; \
	} \
	else \
	{ \
	    r = ar / ai ; \
	    (s) = ai * sqrt (1.0 + r*r) ; \
	} \
    } \
}
#endif

/* -------------------------------------------------------------------------- */

/* print an entry (avoid printing "-0" for negative zero).  */
#define PRINT_ENTRY(a) \
{ \
    if (SCALAR_IS_NONZERO ((a).Real)) \
    { \
	PRINTF ((" (%g", (a).Real)) ; \
    } \
    else \
    { \
	PRINTF ((" (0")) ; \
    } \
    if (SCALAR_IS_LTZERO ((a).Imag)) \
    { \
	PRINTF ((" - %gi)", -(a).Imag)) ; \
    } \
    else if (SCALAR_IS_ZERO ((a).Imag)) \
    { \
	PRINTF ((" + 0i)")) ; \
    } \
    else \
    { \
	PRINTF ((" + %gi)", (a).Imag)) ; \
    } \
}

/* -------------------------------------------------------------------------- */

#endif	/* #ifndef COMPLEX */

/* -------------------------------------------------------------------------- */
/* Double precision, with int's as integers */
/* -------------------------------------------------------------------------- */

#ifdef DINT

#define UMF_analyze		 umfdi_analyze
#define UMF_apply_order		 umfdi_apply_order
#define UMF_assemble		 umfdi_assemble
#define UMF_blas3_update	 umfdi_blas3_update
#define UMF_build_tuples	 umfdi_build_tuples
#define UMF_build_tuples_usage	 umfdi_build_tuples_usage
#define UMF_colamd		 umfdi_colamd
#define UMF_colamd_set_defaults	 umfdi_colamd_set_defaults
#define UMF_create_element	 umfdi_create_element
#define UMF_extend_front	 umfdi_extend_front
#define UMF_free		 umfdi_free
#define UMF_garbage_collection	 umfdi_garbage_collection
#define UMF_get_memory		 umfdi_get_memory
#define UMF_init_front		 umfdi_init_front
#define UMF_is_permutation	 umfdi_is_permutation
#define UMF_kernel		 umfdi_kernel
#define UMF_kernel_init		 umfdi_kernel_init
#define UMF_kernel_init_usage	 umfdi_kernel_init_usage
#define UMF_kernel_wrapup	 umfdi_kernel_wrapup
#define UMF_local_search	 umfdi_local_search
#define UMF_lsolve		 umfdi_lsolve
#define UMF_ltsolve		 umfdi_ltsolve
#define UMF_lhsolve		 umfdi_lhsolve
#define UMF_malloc		 umfdi_malloc
#define UMF_mem_alloc_element	 umfdi_mem_alloc_element
#define UMF_mem_alloc_head_block umfdi_mem_alloc_head_block
#define UMF_mem_alloc_tail_block umfdi_mem_alloc_tail_block
#define UMF_mem_free_tail_block	 umfdi_mem_free_tail_block
#define UMF_mem_init_memoryspace umfdi_mem_init_memoryspace
#define UMF_order_front_tree	 umfdi_order_front_tree
#define UMF_realloc		 umfdi_realloc
#define UMF_report_perm		 umfdi_report_perm
#define UMF_report_vector	 umfdi_report_vector
#define UMF_row_search		 umfdi_row_search
#define UMF_scale_column	 umfdi_scale_column
#define UMF_set_stats		 umfdi_set_stats
#define UMF_solve		 umfdi_solve
#define UMF_symbolic_usage	 umfdi_symbolic_usage
#define UMF_transpose		 umfdi_transpose
#define UMF_tuple_lengths	 umfdi_tuple_lengths
#define UMF_usolve		 umfdi_usolve
#define UMF_utsolve		 umfdi_utsolve
#define UMF_uhsolve		 umfdi_uhsolve
#define UMF_valid_numeric	 umfdi_valid_numeric
#define UMF_valid_symbolic	 umfdi_valid_symbolic
#define UMF_triplet_map_x	 umfdi_triplet_map_x
#define UMF_triplet_map_nox	 umfdi_triplet_map_nox
#define UMF_triplet_nomap_x	 umfdi_triplet_nomap_x
#define UMF_triplet_nomap_nox	 umfdi_triplet_nomap_nox
#define UMFPACK_col_to_triplet	 umfpack_di_col_to_triplet
#define UMFPACK_defaults	 umfpack_di_defaults
#define UMFPACK_free_numeric	 umfpack_di_free_numeric
#define UMFPACK_free_symbolic	 umfpack_di_free_symbolic
#define UMFPACK_get_lunz	 umfpack_di_get_lunz
#define UMFPACK_get_numeric	 umfpack_di_get_numeric
#define UMFPACK_get_symbolic	 umfpack_di_get_symbolic
#define UMFPACK_numeric		 umfpack_di_numeric
#define UMFPACK_qsymbolic	 umfpack_di_qsymbolic
#define UMFPACK_report_control	 umfpack_di_report_control
#define UMFPACK_report_info	 umfpack_di_report_info
#define UMFPACK_report_matrix	 umfpack_di_report_matrix
#define UMFPACK_report_numeric	 umfpack_di_report_numeric
#define UMFPACK_report_perm	 umfpack_di_report_perm
#define UMFPACK_report_status	 umfpack_di_report_status
#define UMFPACK_report_symbolic	 umfpack_di_report_symbolic
#define UMFPACK_report_triplet	 umfpack_di_report_triplet
#define UMFPACK_report_vector	 umfpack_di_report_vector
#define UMFPACK_solve		 umfpack_di_solve
#define UMFPACK_symbolic	 umfpack_di_symbolic
#define UMFPACK_transpose	 umfpack_di_transpose
#define UMFPACK_triplet_to_col	 umfpack_di_triplet_to_col
#define UMFPACK_wsolve		 umfpack_di_wsolve

/* for debugging only: */
#define UMF_malloc_count	 umfdi_malloc_count
#define UMF_debug		 umfdi_debug
#define UMF_nbug		 umfdi_nbug
#define UMF_fbug		 umfdi_fbug
#define UMF_allocfail		 umfdi_allocfail
#define UMF_gprob		 umfdi_gprob
#define UMF_dump_dense		 umfdi_dump_dense
#define UMF_dump_element	 umfdi_dump_element
#define UMF_dump_rowcol		 umfdi_dump_rowcol
#define UMF_dump_matrix		 umfdi_dump_matrix
#define UMF_dump_current_front	 umfdi_dump_current_front
#define UMF_dump_lu		 umfdi_dump_lu
#define UMF_dump_memory		 umfdi_dump_memory
#define UMF_dump_packed_memory	 umfdi_dump_packed_memory
#define UMF_dump_col_matrix	 umfdi_dump_col_matrix
#define UMF_dump_chain		 umfdi_dump_chain
#define UMF_dump_start		 umfdi_dump_start
#define UMF_dump_rowmerge	 umfdi_dump_rowmerge

#endif


/* -------------------------------------------------------------------------- */
/* Double precision, with long's as integers */
/* -------------------------------------------------------------------------- */

#ifdef DLONG

#define UMF_analyze		 umfdl_analyze
#define UMF_apply_order		 umfdl_apply_order
#define UMF_assemble		 umfdl_assemble
#define UMF_blas3_update	 umfdl_blas3_update
#define UMF_build_tuples	 umfdl_build_tuples
#define UMF_build_tuples_usage	 umfdl_build_tuples_usage
#define UMF_colamd		 umfdl_colamd
#define UMF_colamd_set_defaults	 umfdl_colamd_set_defaults
#define UMF_create_element	 umfdl_create_element
#define UMF_extend_front	 umfdl_extend_front
#define UMF_free		 umfdl_free
#define UMF_garbage_collection	 umfdl_garbage_collection
#define UMF_get_memory		 umfdl_get_memory
#define UMF_init_front		 umfdl_init_front
#define UMF_is_permutation	 umfdl_is_permutation
#define UMF_kernel		 umfdl_kernel
#define UMF_kernel_init		 umfdl_kernel_init
#define UMF_kernel_init_usage	 umfdl_kernel_init_usage
#define UMF_kernel_wrapup	 umfdl_kernel_wrapup
#define UMF_local_search	 umfdl_local_search
#define UMF_lsolve		 umfdl_lsolve
#define UMF_ltsolve		 umfdl_ltsolve
#define UMF_lhsolve		 umfdl_lhsolve
#define UMF_malloc		 umfdl_malloc
#define UMF_mem_alloc_element	 umfdl_mem_alloc_element
#define UMF_mem_alloc_head_block umfdl_mem_alloc_head_block
#define UMF_mem_alloc_tail_block umfdl_mem_alloc_tail_block
#define UMF_mem_free_tail_block	 umfdl_mem_free_tail_block
#define UMF_mem_init_memoryspace umfdl_mem_init_memoryspace
#define UMF_order_front_tree	 umfdl_order_front_tree
#define UMF_realloc		 umfdl_realloc
#define UMF_report_perm		 umfdl_report_perm
#define UMF_report_vector	 umfdl_report_vector
#define UMF_row_search		 umfdl_row_search
#define UMF_scale_column	 umfdl_scale_column
#define UMF_set_stats		 umfdl_set_stats
#define UMF_solve		 umfdl_solve
#define UMF_symbolic_usage	 umfdl_symbolic_usage
#define UMF_transpose		 umfdl_transpose
#define UMF_tuple_lengths	 umfdl_tuple_lengths
#define UMF_usolve		 umfdl_usolve
#define UMF_utsolve		 umfdl_utsolve
#define UMF_uhsolve		 umfdl_uhsolve
#define UMF_valid_numeric	 umfdl_valid_numeric
#define UMF_valid_symbolic	 umfdl_valid_symbolic
#define UMF_triplet_map_x	 umfdl_triplet_map_x
#define UMF_triplet_map_nox	 umfdl_triplet_map_nox
#define UMF_triplet_nomap_x	 umfdl_triplet_nomap_x
#define UMF_triplet_nomap_nox	 umfdl_triplet_nomap_nox
#define UMFPACK_col_to_triplet	 umfpack_dl_col_to_triplet
#define UMFPACK_defaults	 umfpack_dl_defaults
#define UMFPACK_free_numeric	 umfpack_dl_free_numeric
#define UMFPACK_free_symbolic	 umfpack_dl_free_symbolic
#define UMFPACK_get_lunz	 umfpack_dl_get_lunz
#define UMFPACK_get_numeric	 umfpack_dl_get_numeric
#define UMFPACK_get_symbolic	 umfpack_dl_get_symbolic
#define UMFPACK_numeric		 umfpack_dl_numeric
#define UMFPACK_qsymbolic	 umfpack_dl_qsymbolic
#define UMFPACK_report_control	 umfpack_dl_report_control
#define UMFPACK_report_info	 umfpack_dl_report_info
#define UMFPACK_report_matrix	 umfpack_dl_report_matrix
#define UMFPACK_report_numeric	 umfpack_dl_report_numeric
#define UMFPACK_report_perm	 umfpack_dl_report_perm
#define UMFPACK_report_status	 umfpack_dl_report_status
#define UMFPACK_report_symbolic	 umfpack_dl_report_symbolic
#define UMFPACK_report_triplet	 umfpack_dl_report_triplet
#define UMFPACK_report_vector	 umfpack_dl_report_vector
#define UMFPACK_solve		 umfpack_dl_solve
#define UMFPACK_symbolic	 umfpack_dl_symbolic
#define UMFPACK_transpose	 umfpack_dl_transpose
#define UMFPACK_triplet_to_col	 umfpack_dl_triplet_to_col
#define UMFPACK_wsolve		 umfpack_dl_wsolve

/* for debugging only: */
#define UMF_malloc_count	 umfdl_malloc_count
#define UMF_debug		 umfdl_debug
#define UMF_nbug		 umfdl_nbug
#define UMF_fbug		 umfdl_fbug
#define UMF_allocfail		 umfdl_allocfail
#define UMF_gprob		 umfdl_gprob
#define UMF_dump_dense		 umfdl_dump_dense
#define UMF_dump_element	 umfdl_dump_element
#define UMF_dump_rowcol		 umfdl_dump_rowcol
#define UMF_dump_matrix		 umfdl_dump_matrix
#define UMF_dump_current_front	 umfdl_dump_current_front
#define UMF_dump_lu		 umfdl_dump_lu
#define UMF_dump_memory		 umfdl_dump_memory
#define UMF_dump_packed_memory	 umfdl_dump_packed_memory
#define UMF_dump_col_matrix	 umfdl_dump_col_matrix
#define UMF_dump_chain		 umfdl_dump_chain
#define UMF_dump_start		 umfdl_dump_start
#define UMF_dump_rowmerge	 umfdl_dump_rowmerge

#endif

/* -------------------------------------------------------------------------- */
/* Complex double precision, with int's as integers */
/* -------------------------------------------------------------------------- */

#ifdef ZINT

#define UMF_analyze		 umfzi_analyze
#define UMF_apply_order		 umfzi_apply_order
#define UMF_assemble		 umfzi_assemble
#define UMF_blas3_update	 umfzi_blas3_update
#define UMF_build_tuples	 umfzi_build_tuples
#define UMF_build_tuples_usage	 umfzi_build_tuples_usage
#define UMF_colamd		 umfzi_colamd
#define UMF_colamd_set_defaults	 umfzi_colamd_set_defaults
#define UMF_create_element	 umfzi_create_element
#define UMF_extend_front	 umfzi_extend_front
#define UMF_free		 umfzi_free
#define UMF_garbage_collection	 umfzi_garbage_collection
#define UMF_get_memory		 umfzi_get_memory
#define UMF_init_front		 umfzi_init_front
#define UMF_is_permutation	 umfzi_is_permutation
#define UMF_kernel		 umfzi_kernel
#define UMF_kernel_init		 umfzi_kernel_init
#define UMF_kernel_init_usage	 umfzi_kernel_init_usage
#define UMF_kernel_wrapup	 umfzi_kernel_wrapup
#define UMF_local_search	 umfzi_local_search
#define UMF_lsolve		 umfzi_lsolve
#define UMF_ltsolve		 umfzi_ltsolve
#define UMF_lhsolve		 umfzi_lhsolve
#define UMF_malloc		 umfzi_malloc
#define UMF_mem_alloc_element	 umfzi_mem_alloc_element
#define UMF_mem_alloc_head_block umfzi_mem_alloc_head_block
#define UMF_mem_alloc_tail_block umfzi_mem_alloc_tail_block
#define UMF_mem_free_tail_block	 umfzi_mem_free_tail_block
#define UMF_mem_init_memoryspace umfzi_mem_init_memoryspace
#define UMF_order_front_tree	 umfzi_order_front_tree
#define UMF_realloc		 umfzi_realloc
#define UMF_report_perm		 umfzi_report_perm
#define UMF_report_vector	 umfzi_report_vector
#define UMF_row_search		 umfzi_row_search
#define UMF_scale_column	 umfzi_scale_column
#define UMF_set_stats		 umfzi_set_stats
#define UMF_solve		 umfzi_solve
#define UMF_symbolic_usage	 umfzi_symbolic_usage
#define UMF_transpose		 umfzi_transpose
#define UMF_tuple_lengths	 umfzi_tuple_lengths
#define UMF_usolve		 umfzi_usolve
#define UMF_utsolve		 umfzi_utsolve
#define UMF_uhsolve		 umfzi_uhsolve
#define UMF_valid_numeric	 umfzi_valid_numeric
#define UMF_valid_symbolic	 umfzi_valid_symbolic
#define UMF_triplet_map_x	 umfzi_triplet_map_x
#define UMF_triplet_map_nox	 umfzi_triplet_map_nox
#define UMF_triplet_nomap_x	 umfzi_triplet_nomap_x
#define UMF_triplet_nomap_nox	 umfzi_triplet_nomap_nox
#define UMFPACK_col_to_triplet	 umfpack_zi_col_to_triplet
#define UMFPACK_defaults	 umfpack_zi_defaults
#define UMFPACK_free_numeric	 umfpack_zi_free_numeric
#define UMFPACK_free_symbolic	 umfpack_zi_free_symbolic
#define UMFPACK_get_lunz	 umfpack_zi_get_lunz
#define UMFPACK_get_numeric	 umfpack_zi_get_numeric
#define UMFPACK_get_symbolic	 umfpack_zi_get_symbolic
#define UMFPACK_numeric		 umfpack_zi_numeric
#define UMFPACK_qsymbolic	 umfpack_zi_qsymbolic
#define UMFPACK_report_control	 umfpack_zi_report_control
#define UMFPACK_report_info	 umfpack_zi_report_info
#define UMFPACK_report_matrix	 umfpack_zi_report_matrix
#define UMFPACK_report_numeric	 umfpack_zi_report_numeric
#define UMFPACK_report_perm	 umfpack_zi_report_perm
#define UMFPACK_report_status	 umfpack_zi_report_status
#define UMFPACK_report_symbolic	 umfpack_zi_report_symbolic
#define UMFPACK_report_triplet	 umfpack_zi_report_triplet
#define UMFPACK_report_vector	 umfpack_zi_report_vector
#define UMFPACK_solve		 umfpack_zi_solve
#define UMFPACK_symbolic	 umfpack_zi_symbolic
#define UMFPACK_transpose	 umfpack_zi_transpose
#define UMFPACK_triplet_to_col	 umfpack_zi_triplet_to_col
#define UMFPACK_wsolve		 umfpack_zi_wsolve

/* for debugging only: */
#define UMF_malloc_count	 umfzi_malloc_count
#define UMF_debug		 umfzi_debug
#define UMF_nbug		 umfzi_nbug
#define UMF_fbug		 umfzi_fbug
#define UMF_allocfail		 umfzi_allocfail
#define UMF_gprob		 umfzi_gprob
#define UMF_dump_dense		 umfzi_dump_dense
#define UMF_dump_element	 umfzi_dump_element
#define UMF_dump_rowcol		 umfzi_dump_rowcol
#define UMF_dump_matrix		 umfzi_dump_matrix
#define UMF_dump_current_front	 umfzi_dump_current_front
#define UMF_dump_lu		 umfzi_dump_lu
#define UMF_dump_memory		 umfzi_dump_memory
#define UMF_dump_packed_memory	 umfzi_dump_packed_memory
#define UMF_dump_col_matrix	 umfzi_dump_col_matrix
#define UMF_dump_chain		 umfzi_dump_chain
#define UMF_dump_start		 umfzi_dump_start
#define UMF_dump_rowmerge	 umfzi_dump_rowmerge

#endif


/* -------------------------------------------------------------------------- */
/* Complex double precision, with long's as integers */
/* -------------------------------------------------------------------------- */

#ifdef ZLONG

#define UMF_analyze		 umfzl_analyze
#define UMF_apply_order		 umfzl_apply_order
#define UMF_assemble		 umfzl_assemble
#define UMF_blas3_update	 umfzl_blas3_update
#define UMF_build_tuples	 umfzl_build_tuples
#define UMF_build_tuples_usage	 umfzl_build_tuples_usage
#define UMF_colamd		 umfzl_colamd
#define UMF_colamd_set_defaults	 umfzl_colamd_set_defaults
#define UMF_create_element	 umfzl_create_element
#define UMF_extend_front	 umfzl_extend_front
#define UMF_free		 umfzl_free
#define UMF_garbage_collection	 umfzl_garbage_collection
#define UMF_get_memory		 umfzl_get_memory
#define UMF_init_front		 umfzl_init_front
#define UMF_is_permutation	 umfzl_is_permutation
#define UMF_kernel		 umfzl_kernel
#define UMF_kernel_init		 umfzl_kernel_init
#define UMF_kernel_init_usage	 umfzl_kernel_init_usage
#define UMF_kernel_wrapup	 umfzl_kernel_wrapup
#define UMF_local_search	 umfzl_local_search
#define UMF_lsolve		 umfzl_lsolve
#define UMF_ltsolve		 umfzl_ltsolve
#define UMF_lhsolve		 umfzl_lhsolve
#define UMF_malloc		 umfzl_malloc
#define UMF_mem_alloc_element	 umfzl_mem_alloc_element
#define UMF_mem_alloc_head_block umfzl_mem_alloc_head_block
#define UMF_mem_alloc_tail_block umfzl_mem_alloc_tail_block
#define UMF_mem_free_tail_block	 umfzl_mem_free_tail_block
#define UMF_mem_init_memoryspace umfzl_mem_init_memoryspace
#define UMF_order_front_tree	 umfzl_order_front_tree
#define UMF_realloc		 umfzl_realloc
#define UMF_report_perm		 umfzl_report_perm
#define UMF_report_vector	 umfzl_report_vector
#define UMF_row_search		 umfzl_row_search
#define UMF_scale_column	 umfzl_scale_column
#define UMF_set_stats		 umfzl_set_stats
#define UMF_solve		 umfzl_solve
#define UMF_symbolic_usage	 umfzl_symbolic_usage
#define UMF_transpose		 umfzl_transpose
#define UMF_tuple_lengths	 umfzl_tuple_lengths
#define UMF_usolve		 umfzl_usolve
#define UMF_utsolve		 umfzl_utsolve
#define UMF_uhsolve		 umfzl_uhsolve
#define UMF_valid_numeric	 umfzl_valid_numeric
#define UMF_valid_symbolic	 umfzl_valid_symbolic
#define UMF_triplet_map_x	 umfzl_triplet_map_x
#define UMF_triplet_map_nox	 umfzl_triplet_map_nox
#define UMF_triplet_nomap_x	 umfzl_triplet_nomap_x
#define UMF_triplet_nomap_nox	 umfzl_triplet_nomap_nox
#define UMFPACK_col_to_triplet	 umfpack_zl_col_to_triplet
#define UMFPACK_defaults	 umfpack_zl_defaults
#define UMFPACK_free_numeric	 umfpack_zl_free_numeric
#define UMFPACK_free_symbolic	 umfpack_zl_free_symbolic
#define UMFPACK_get_lunz	 umfpack_zl_get_lunz
#define UMFPACK_get_numeric	 umfpack_zl_get_numeric
#define UMFPACK_get_symbolic	 umfpack_zl_get_symbolic
#define UMFPACK_numeric		 umfpack_zl_numeric
#define UMFPACK_qsymbolic	 umfpack_zl_qsymbolic
#define UMFPACK_report_control	 umfpack_zl_report_control
#define UMFPACK_report_info	 umfpack_zl_report_info
#define UMFPACK_report_matrix	 umfpack_zl_report_matrix
#define UMFPACK_report_numeric	 umfpack_zl_report_numeric
#define UMFPACK_report_perm	 umfpack_zl_report_perm
#define UMFPACK_report_status	 umfpack_zl_report_status
#define UMFPACK_report_symbolic	 umfpack_zl_report_symbolic
#define UMFPACK_report_triplet	 umfpack_zl_report_triplet
#define UMFPACK_report_vector	 umfpack_zl_report_vector
#define UMFPACK_solve		 umfpack_zl_solve
#define UMFPACK_symbolic	 umfpack_zl_symbolic
#define UMFPACK_transpose	 umfpack_zl_transpose
#define UMFPACK_triplet_to_col	 umfpack_zl_triplet_to_col
#define UMFPACK_wsolve		 umfpack_zl_wsolve

/* for debugging only: */
#define UMF_malloc_count	 umfzl_malloc_count
#define UMF_debug		 umfzl_debug
#define UMF_nbug		 umfzl_nbug
#define UMF_fbug		 umfzl_fbug
#define UMF_allocfail		 umfzl_allocfail
#define UMF_gprob		 umfzl_gprob
#define UMF_dump_dense		 umfzl_dump_dense
#define UMF_dump_element	 umfzl_dump_element
#define UMF_dump_rowcol		 umfzl_dump_rowcol
#define UMF_dump_matrix		 umfzl_dump_matrix
#define UMF_dump_current_front	 umfzl_dump_current_front
#define UMF_dump_lu		 umfzl_dump_lu
#define UMF_dump_memory		 umfzl_dump_memory
#define UMF_dump_packed_memory	 umfzl_dump_packed_memory
#define UMF_dump_col_matrix	 umfzl_dump_col_matrix
#define UMF_dump_chain		 umfzl_dump_chain
#define UMF_dump_start		 umfzl_dump_start
#define UMF_dump_rowmerge	 umfzl_dump_rowmerge

#endif
