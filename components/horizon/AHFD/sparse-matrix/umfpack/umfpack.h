/* ========================================================================== */
/* === umfpack.h ============================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    This is the umfpack.h include file, and should be included in all user code
    that uses UMFPACK.  Do not include any of the umf_* header files in user
    code.  All routines in UMFPACK starting with "umfpack_" are user-callable.
    All other routines are prefixed "umfXY_", (where X is d or z, and Y is
    i or l) and are not user-callable.
*/

#ifndef UMFPACK_H
#define UMFPACK_H

/* -------------------------------------------------------------------------- */
/* size of Info and Control arrays */
/* -------------------------------------------------------------------------- */

#define UMFPACK_INFO 90
#define UMFPACK_CONTROL 20

/* -------------------------------------------------------------------------- */
/* User-callable routines */
/* -------------------------------------------------------------------------- */

/* Primary routines: */
#include "umfpack_symbolic.h"
#include "umfpack_numeric.h"
#include "umfpack_solve.h"
#include "umfpack_free_symbolic.h"
#include "umfpack_free_numeric.h"

/* Alternative routines: */
#include "umfpack_defaults.h"
#include "umfpack_qsymbolic.h"
#include "umfpack_wsolve.h"

/* Matrix manipulation routines: */
#include "umfpack_triplet_to_col.h"
#include "umfpack_col_to_triplet.h"
#include "umfpack_transpose.h"

/* Getting the contents of the Symbolic and Numeric opaque objects: */
#include "umfpack_get_lunz.h"
#include "umfpack_get_numeric.h"
#include "umfpack_get_symbolic.h"

/* Reporting routines (the above 14 routines print nothing): */
#include "umfpack_report_status.h"
#include "umfpack_report_info.h"
#include "umfpack_report_control.h"
#include "umfpack_report_matrix.h"
#include "umfpack_report_triplet.h"
#include "umfpack_report_vector.h"
#include "umfpack_report_symbolic.h"
#include "umfpack_report_numeric.h"
#include "umfpack_report_perm.h"

/* Utility routines: */
#include "umfpack_timer.h"


/* -------------------------------------------------------------------------- */
/* Version, copyright, and license */
/* -------------------------------------------------------------------------- */

#define UMFPACK_VERSION "UMFPACK V4.0 (Apr 11, 2002)"

#define UMFPACK_COPYRIGHT \
"UMFPACK:  Copyright (c) 2002 by Timothy A. Davis.  All Rights Reserved.\n"

#define UMFPACK_LICENSE_PART1 \
"\nUMFPACK License:\n" \
"\n" \
"   Your use or distribution of UMFPACK or any modified version of\n" \
"   UMFPACK implies that you agree to this License.\n" \
"\n" \
"   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY\n" \
"   EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.\n"
#define UMFPACK_LICENSE_PART2 \
"\n" \
"   Permission is hereby granted to use or copy this program, provided\n" \
"   that the Copyright, this License, and the Availability of the original\n" \
"   version is retained on all copies.  User documentation of any code that\n" \
"   uses UMFPACK or any modified version of UMFPACK code must cite the\n" \
"   Copyright, this License, the Availability note, and \"Used by permission.\"\n"
#define UMFPACK_LICENSE_PART3 \
"   Permission to modify the code and to distribute modified code is granted,\n" \
"   provided the Copyright, this License, and the Availability note are\n" \
"   retained, and a notice that the code was modified is included.  This\n" \
"   software was developed with support from the National Science Foundation,\n" \
"   and is provided to you free of charge.\n" \
"\n" \
"Availability: http://www.cise.ufl.edu/research/sparse/umfpack\n" \
"\n"

/* -------------------------------------------------------------------------- */
/* contents of Info */
/* -------------------------------------------------------------------------- */

/* Note that umfpack_report.m must coincide with these definitions. */

/* returned by all routines that use Info: */
#define UMFPACK_STATUS 0
#define UMFPACK_NROW 1
#define UMFPACK_NCOL 16
#define UMFPACK_NZ 2

/* computed in UMFPACK_*symbolic and UMFPACK_numeric: */
#define UMFPACK_SIZE_OF_UNIT 3

/* computed in UMFPACK_*symbolic: */
#define UMFPACK_SIZE_OF_INT 4
#define UMFPACK_SIZE_OF_LONG 5
#define UMFPACK_SIZE_OF_POINTER 6
#define UMFPACK_SIZE_OF_ENTRY 7
#define UMFPACK_NDENSE_ROW 8
#define UMFPACK_NEMPTY_ROW 9
#define UMFPACK_NDENSE_COL 10
#define UMFPACK_NEMPTY_COL 11
#define UMFPACK_SYMBOLIC_DEFRAG 12
#define UMFPACK_SYMBOLIC_PEAK_MEMORY 13
#define UMFPACK_SYMBOLIC_SIZE 14
#define UMFPACK_SYMBOLIC_TIME 15

/* Info [17..19] unused */

/* estimates computed in UMFPACK_*symbolic: */
#define UMFPACK_NUMERIC_SIZE_ESTIMATE 20
#define UMFPACK_PEAK_MEMORY_ESTIMATE 21
#define UMFPACK_FLOPS_ESTIMATE 22
#define UMFPACK_LNZ_ESTIMATE 23
#define UMFPACK_UNZ_ESTIMATE 24
#define UMFPACK_VARIABLE_INIT_ESTIMATE 25
#define UMFPACK_VARIABLE_PEAK_ESTIMATE 26
#define UMFPACK_VARIABLE_FINAL_ESTIMATE 27
#define UMFPACK_MAX_FRONT_SIZE_ESTIMATE 28

/* Info [29..39] unused */

/* exact values, (estimates shown above) computed in UMFPACK_numeric: */
#define UMFPACK_NUMERIC_SIZE 40
#define UMFPACK_PEAK_MEMORY 41
#define UMFPACK_FLOPS 42
#define UMFPACK_LNZ 43
#define UMFPACK_UNZ 44
#define UMFPACK_VARIABLE_INIT 45
#define UMFPACK_VARIABLE_PEAK 46
#define UMFPACK_VARIABLE_FINAL 47
#define UMFPACK_MAX_FRONT_SIZE 48

/* Info [49..59] unused */

/* computed in UMFPACK_numeric: */
#define UMFPACK_NUMERIC_DEFRAG 60
#define UMFPACK_NUMERIC_REALLOC 61
#define UMFPACK_NUMERIC_COSTLY_REALLOC 62
#define UMFPACK_COMPRESSED_PATTERN 63
#define UMFPACK_LU_ENTRIES 64
#define UMFPACK_NUMERIC_TIME 65
#define UMFPACK_UDIAG_NZ 66
#define UMFPACK_RCOND 67

/* Info [68..79] unused */

/* computed in UMFPACK_solve: */
#define UMFPACK_IR_TAKEN 80
#define UMFPACK_IR_ATTEMPTED 81
#define UMFPACK_OMEGA1 82
#define UMFPACK_OMEGA2 83
#define UMFPACK_SOLVE_FLOPS 84
#define UMFPACK_SOLVE_TIME 85

/* Info [86..89] unused */

/* Unused parts of Info may be used in future versions of UMFPACK. */


/* -------------------------------------------------------------------------- */
/* contents of Control */
/* -------------------------------------------------------------------------- */

/* used in all UMFPACK_report_* routines: */
#define UMFPACK_PRL 0

/* used in UMFPACK_*symbolic only: */
#define UMFPACK_DENSE_ROW 1
#define UMFPACK_DENSE_COL 2

/* used in UMFPACK_numeric only: */
#define UMFPACK_PIVOT_TOLERANCE 3
#define UMFPACK_BLOCK_SIZE 4
#define UMFPACK_RELAXED_AMALGAMATION 5
#define UMFPACK_ALLOC_INIT 6
/* #define UMFPACK_PIVOT_OPTION 12: obsolete */
#define UMFPACK_RELAXED2_AMALGAMATION 13
#define UMFPACK_RELAXED3_AMALGAMATION 14

/* used in UMFPACK_*solve only: */
#define UMFPACK_IRSTEP 7

/* compile-time settings - Control [8..11] cannot be changed at run time: */
#define UMFPACK_COMPILED_WITH_BLAS 8
#define UMFPACK_COMPILED_FOR_MATLAB 9
#define UMFPACK_COMPILED_WITH_GETRUSAGE 10
#define UMFPACK_COMPILED_IN_DEBUG_MODE 11

/* Control [12, 15...19] unused */

/* Unused parts of Control may be used in future versions of UMFPACK. */


/* -------------------------------------------------------------------------- */
/* default values of Control [0..7,13..14]: */
/* -------------------------------------------------------------------------- */

/* Note that the default block sized changed for Version 3.1 and following. */

#define UMFPACK_DEFAULT_PRL 1
#define UMFPACK_DEFAULT_DENSE_ROW 0.2
#define UMFPACK_DEFAULT_DENSE_COL 0.2
#define UMFPACK_DEFAULT_PIVOT_TOLERANCE 0.1
#define UMFPACK_DEFAULT_BLOCK_SIZE 24
#define UMFPACK_DEFAULT_RELAXED_AMALGAMATION 0.25
#define UMFPACK_DEFAULT_RELAXED2_AMALGAMATION 0.1
#define UMFPACK_DEFAULT_RELAXED3_AMALGAMATION 0.125
#define UMFPACK_DEFAULT_ALLOC_INIT 0.7
#define UMFPACK_DEFAULT_IRSTEP 2
/* #define UMFPACK_DEFAULT_PIVOT_OPTION 0: obsolete */

/* default values of Control [0..7,13..14] may change in future versions */
/* of UMFPACK. */

/* -------------------------------------------------------------------------- */
/* status codes */
/* -------------------------------------------------------------------------- */

#define UMFPACK_OK (0)

/* status > 0 means a warning, but the method was successful anyway. */
/* A Symbolic or Numeric object was still created. */
#define UMFPACK_WARNING_singular_matrix (1)

/* status < 0 means an error, and the method was not successful. */
/* No Symbolic of Numeric object was created. */
#define UMFPACK_ERROR_out_of_memory (-1)
#define UMFPACK_ERROR_invalid_Numeric_object (-3)
#define UMFPACK_ERROR_invalid_Symbolic_object (-4)
#define UMFPACK_ERROR_argument_missing (-5)
#define UMFPACK_ERROR_n_nonpositive (-6)
#define UMFPACK_ERROR_nz_negative (-7)
#define UMFPACK_ERROR_jumbled_matrix (-8)
#define UMFPACK_ERROR_Ap0_nonzero (-9)
#define UMFPACK_ERROR_row_index_out_of_bounds (-10)
#define UMFPACK_ERROR_different_pattern (-11)
#define UMFPACK_ERROR_col_length_negative (-12)
#define UMFPACK_ERROR_invalid_system (-13)
#define UMFPACK_ERROR_invalid_triplet (-14)
#define UMFPACK_ERROR_invalid_permutation (-15)
#define UMFPACK_ERROR_problem_too_large (-16)
#define UMFPACK_ERROR_internal_error (-911)

/* -------------------------------------------------------------------------- */
/* solve codes */
/* -------------------------------------------------------------------------- */

/* Solve the system ( )x=b, where ( ) is defined below.  "t" refers to the */
/* linear algebraic transpose (complex conjugate if A is complex), or the (') */
/* operator in MATLAB.  "at" refers to the array transpose, or the (.') */
/* operator in MATLAB. */

#define UMFPACK_A	(0)	/* Ax=b		*/
#define UMFPACK_At	(1)	/* A'x=b	*/
#define UMFPACK_Aat	(2)	/* A.'x=b	*/

#define UMFPACK_Pt_L	(3)	/* P'Lx=b	*/
#define UMFPACK_L	(4)	/* Lx=b		*/
#define UMFPACK_Lt_P	(5)	/* L'Px=b	*/
#define UMFPACK_Lat_P	(6)	/* L.'Px=b	*/
#define UMFPACK_Lt	(7)	/* L'x=b	*/
#define UMFPACK_Lat	(8)	/* L.'x=b	*/

#define UMFPACK_U_Qt	(9)	/* UQ'x=b	*/
#define UMFPACK_U	(10)	/* Ux=b		*/
#define UMFPACK_Q_Ut	(11)	/* QU'x=b	*/
#define UMFPACK_Q_Uat	(12)	/* QU.'x=b	*/
#define UMFPACK_Ut	(13)	/* U'x=b	*/
#define UMFPACK_Uat	(14)	/* U.'x=b	*/

/* -------------------------------------------------------------------------- */

/* Integer constants are used for status and solve codes instead of enum */
/* to make it easier for a Fortran code to call UMFPACK. */

#endif /* UMFPACK_H */
