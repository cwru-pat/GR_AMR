/* ========================================================================== */
/* === umf_config.h ========================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    This file controls the compile-time configuration of UMFPACK.  Modify the
    Makefile, the architecture-dependent Make.* file, and this file if
    necessary, to control these options.  The following flags may be given
    as options to your C compiler (as in "cc -DNBLAS", for example).  These
    flags are normally placed in your CONFIG string, defined in your Make.*.

	-DNBLAS

	    BLAS mode.  If -DNBLAS is set, then no BLAS will be used.  Vanilla
	    C code will be used instead.  This is portable, and easier to
	    install, but you won't get the best performance.

	    If -DNBLAS is not set, then externally-available BLAS routines
	    (dgemm, dger, and dgemv or the equivalent C-BLAS routines) will be
	    used.  This will give you the best performance, but perhaps at the
	    expense of portability.

	    The default is to use the BLAS, for both the C-callable umfpack.a
	    library and the MATLAB mexFunction.  If you have trouble installing
	    UMFPACK, set -DNBLAS.

	-DNCBLAS

	    If -DNCBLAS is set, then the C-BLAS will not be called.  This is the
	    default when compiling the MATLAB mexFunction, or when compiling
	    umfpack.a on Sun Solaris or SGI IRIX.

	    If -DNCBLAS is not set, then the C-BLAS interface to the BLAS is
	    used.  If your vendor-supplied BLAS library does not have a C-BLAS
	    interface, you can obtain the ATLAS BLAS, available at
	    http://www.netlib.org/atlas.

	    Using the C-BLAS is the default when compiling umfpack.a on all
	    architectures except Sun Solaris (the Sun Performance Library is
	    somewhat faster).  The ANSI C interface to the BLAS is fully
	    portable.

	    This flag is ignored if -DNBLAS is set.

	-DLP64

	    This should be defined if you are compiling in the LP64 model
	    (32 bit int's, 64 bit long's, and 64 bit pointers).  In Solaris,
	    this is obtained with the flags -xtarget=ultra -xarch=v9 for
	    the cc compiler (for example).

	-DLONGBLAS

	    If not defined, then the BLAS are not called in the long integer
	    version of UMFPACK (the umfpack_*l_* routines).  The most common
	    definitions of the BLAS, unfortunately, use int arguments, and
	    are thus not suitable for use in the LP64 model.  Only the Sun
	    Performance Library, as far as I can tell, has a version of the
	    BLAS that allows long integer (64-bit) input arguments.  This
	    flag is set automatically in Sun Solaris if you are using the
	    Sun Performance BLAS.  You can set it yourself, too, if your BLAS
	    routines can take long integer input arguments.

	-DNSUNPERF

	    Applies only to Sun Solaris.  If -DNSUNPERF is set, then the Sun
	    Performance Library BLAS will not be used.

	    The Sun Performance Library BLAS is used by default when compiling
	    the C-callable umfpack.a library on Sun Solaris.

	    This flag is ignored if -DNBLAS is set.

	-DNSCSL

	    Applies only to SGI IRIX.  If -DSCSL is set, then the SGI SCSL
	    Scientific Library BLAS will not be used.

	    The SGI SCSL Scientific Library BLAS is used by default when
	    compiling the C-callable umfpack.a library on SGI IRIX.

	    This flag is ignored if -DNBLAS is set.

	-DGETRUSAGE

	    If -DGETRUSAGE is set, then your system's getrusage routine will be
	    used for getting the process CPU time.  Otherwise the ANSI C clock
	    routine will be used.  The default is to use getrusage on Sun
	    Solaris, SGI Irix, Linux, and AIX (IBM RS 6000) and to use clock on
	    all other architectures.

    You should normally not set these flags yourself:

	-DBLAS_BY_VALUE		if scalars are passed by value, not reference
	-DBLAS_NO_UNDERSCORE	if no underscore should be appended
	-DBLAS_CHAR_ARG		if BLAS options are single char's, not strings

	    The BLAS options are normally set automatically.  If your
	    architecture cannot be determined (see UMFPACK_ARCHITECTURE, below)
	    then you may need to set these flags yourself.

	-DMATLAB_MEX_FILE

	    This flag is turned on when compiling the umfpack mexFunction for
	    use in MATLAB.  When compiling the MATLAB mexFunction, the MATLAB
	    BLAS are used (unless -DNBLAS is set).  The -DNCBLAS, -DNSCSL, and
	    -DNSUNPERF flags are all ignored.

	-DMATHWORKS

	    This flag is turned on when compiling umfpack as a built-in routine
	    in MATLAB.  The MATLAB BLAS are used for all architectures (-DNBLAS,
	    -DNCBLAS, -DNSCSL, and -DNSUNPERF flags are all ignored).  Internal
	    routines utMalloc, utFree, utRealloc, utPrintf, utDivideComplex,
	    and utFdlibm_hypot are used, and the "util.h" file is included.  
	    This avoids the problem discussed in the User Guide regarding memory
	    allocation in MATLAB.  utMalloc returns NULL on failure, instead of
	    terminating the mexFunction (which is what mxMalloc does).  However,
	    the ut* routines are not documented by The MathWorks, Inc., so I
	    cannot guarantee that you will always be able to use them.

	-DNDEBUG

	    Debugging mode (if NDEBUG is not defined).  The default, of course,
	    is no debugging.  Turning on debugging takes some work (see below).
	    If you do not edit this file, then debugging is turned off anyway,
	    regardless of whether or not -DNDEBUG is specified in your compiler
	    options.

    You can change the memory allocator routines by editting the definitions
    of ALLOCATE, FREE, and REALLOCATE, below.
*/



/* ========================================================================== */
/* === NDEBUG =============================================================== */
/* ========================================================================== */

/*
    UMFPACK will be exceedingly slow when running in debug mode.  The next three
    lines ensure that debugging is turned off.  If you want to compile UMFPACK
    in debugging mode, you must comment out the three lines below:
*/
#ifndef NDEBUG
#define NDEBUG
#endif

/*
    Next, you must either remove the -DNDEBUG option in the Makefile, or simply
    add the following line:
#undef NDEBUG
*/

/* ========================================================================== */
/* === Memory allocator ===================================================== */
/* ========================================================================== */

/* The MATLAB mexFunction uses MATLAB's memory manager, while the C-callable */
/* umfpack.a library uses the ANSI C malloc, free, and realloc routines. */

#ifdef MATLAB_MEX_FILE
#define ALLOCATE mxMalloc
#define FREE mxFree
#define REALLOCATE mxRealloc
#else
#ifdef MATHWORKS
/* Compiling UMFPACK as a built-in routine. */
/* Since UMFPACK carefully checks for out-of-memory after every allocation, */
/* we can use ut* routines here. */
#define ALLOCATE utMalloc
#define FREE utFree
#define REALLOCATE utRealloc
#else
#define ALLOCATE malloc
#define FREE free
#define REALLOCATE realloc
#endif
#endif


/* ========================================================================== */
/* === PRINTF macro ========================================================= */
/* ========================================================================== */

/* All output goes through the PRINTF macro.  Printing occurs only from the */
/* UMFPACK_report_* routines. */

#ifdef MATLAB_MEX_FILE
#define PRINTF(params) { (void) mexPrintf params ; }
#else
#ifdef MATHWORKS
/* Already #included "util.h" in umf_version.h */
#define PRINTF(params) { (void) utPrintf params ; }
#else
#define PRINTF(params) { (void) printf params ; }
#endif
#endif


/* ========================================================================== */
/* === 0-based or 1-based printing ========================================== */
/* ========================================================================== */

#if defined (MATLAB_MEX_FILE) && defined (NDEBUG)
/* In MATLAB, matrices are 1-based to the user, but 0-based internally. */
/* One is added to all row and column indices when printing matrices */
/* in UMFPACK_report_*.  */
#define INDEX(i) ((i)+1)
#else
/* In ANSI C, matrices are 0-based and indices are reported as such. */
/* This mode is also used for debug mode, and if MATHWORKS is defined rather */
/* than MATLAB_MEX_FILE. */
#define INDEX(i) (i)
#endif


/* ========================================================================== */
/* === Timer ================================================================ */
/* ========================================================================== */

/*
    If you have the getrusage routine (all Unix systems I've test do), then use
    that.  Otherwise, use the ANSI C clock function.   Note that on many
    systems, the ANSI clock function wraps around after only 2147 seconds, or
    about 36 minutes.  BE CAREFUL:  if you compare the run time of UMFPACK with
    other sparse matrix packages, be sure to use the same timer.  See
    umfpack_timer.c for the timer used by UMFPACK.
*/

/* Sun Solaris, SGI Irix, Linux, Compaq Alpha, and IBM RS 6000 all have */
/* getrusage.  It's in BSD unix, so perhaps all unix systems have it. */
#if defined (UMF_SOL2) || defined (UMF_SGI) || defined (UMF_LINUX) \
|| defined (UMF_ALPHA) || defined (UMF_AIX)
#define GETRUSAGE
#endif


/* ========================================================================== */
/* === BLAS ================================================================= */
/* ========================================================================== */

/*
    The adventure begins.  Figure out how to call the BLAS ...

    This works, but it is incredibly ugly.  The C-BLAS was supposed to solve
    this problem, and make it easier to interface a C program to the BLAS.
    Unfortunately, the C-BLAS does not have a "long" integer (64 bit) version.
    Various vendors have done their own 64-bit BLAS.  Sun has dgemm_64 routines
    with long integers, SGI has a 64-bit dgemm in their scsl_blas_i8 library
    with long long integers, and so on.

    Different vendors also have different ways of defining a complex number,
    some using struct's.  That's a bad idea.  See umf_version.h for the better
    way to do it (the method that was also chosen for the complex C-BLAS).

    To make matters worse, SGI's SCSL BLAS has a C-BLAS interface which
    differs from the ATLAS C-BLAS interface (see immediately below);
    although a more recent version of SGI's C-BLAS interface is correct
    if SCSL_VOID_ARGS is defined.
*/


/* -------------------------------------------------------------------------- */
/* Determine which BLAS to use. */
/* -------------------------------------------------------------------------- */

#if defined (MATHWORKS)
#define USE_MATLAB_BLAS

#elif defined (NBLAS)
#define USE_NO_BLAS

#elif defined (MATLAB_MEX_FILE)
#define USE_MATLAB_BLAS

#elif defined (UMF_SOL2) && !defined (NSUNPERF)
#define USE_SUNPERF_BLAS

#elif defined (UMF_SGI) && !defined (NSCSL)
#define USE_SCSL_BLAS

#elif !defined (NCBLAS)
#define USE_C_BLAS

#else
#define USE_FORTRAN_BLAS
#endif

/* -------------------------------------------------------------------------- */
/* int vs. long integer arguments */
/* -------------------------------------------------------------------------- */

/*
    Determine if the BLAS exists for the long integer version.  It exists if
    LONGBLAS is defined in the Makefile, or if using the BLAS from the
    Sun Performance Library, or SGI's SCSL Scientific Library.
*/

#if defined (USE_SUNPERF_BLAS) || defined (USE_SCSL_BLAS)
#ifndef LONGBLAS
#define LONGBLAS
#endif
#endif

/* do not use the BLAS if Int's are long and LONGBLAS is not defined */
#if defined (LONG_INTEGER) && !defined (LONGBLAS) && !defined (USE_NO_BLAS)
#define USE_NO_BLAS
#endif


/* -------------------------------------------------------------------------- */
/* Use (void *) arguments for the SGI */
/* -------------------------------------------------------------------------- */

#if defined (UMF_SGI)
/*
    Use (void *) pointers for complex types in SCSL.
    The ATLAS C-BLAS, and the SGI C-BLAS differ.  The former uses (void *)
    arguments, the latter uses SCSL_ZOMPLEX_T, which are either scsl_zomplex
    or (void *).  Using (void *) is simpler, and is selected by defining
    SCSL_VOID_ARGS, below.  The cc compiler doesn't complain, but gcc is
    more picky, and generates a warning without this next statement.
    With gcc and the 07/09/98 version of SGI's cblas.h, spurious warnings
    about complex BLAS arguments will be reported anyway.  This is because this
    older version of SGI's cblas.h does not make use of the SCSL_VOID_ARGS
    parameter, which is present in the 12/6/01 version of SGI's cblas.h.  You
    can safely ignore these warnings.
*/
#define SCSL_VOID_ARGS
#endif


/* -------------------------------------------------------------------------- */
/* The BLAS exists, construct appropriate macros */
/* -------------------------------------------------------------------------- */

#if !defined (USE_NO_BLAS)		/* { */

/*
    If the compile-time flag -DNBLAS is defined, then the BLAS are not used,
    portable vanilla C code is used instead, and the remainder of this file
    is ignored.

    Using the BLAS is much faster, but how C calls the Fortran BLAS is
    machine-dependent and thus can cause portability problems.  Thus, use
    -DNBLAS to ensure portability (at the expense of speed).

    Preferences:

	*** The best interface to use, regardless of the option you select
	    below, is the standard C-BLAS interface.  Not all vendor-supplied
	    BLAS libraries use this interface (as of April 2001).  The only
	    problem with this interface is that it does not extend to the LP64
	    model.  The C-BLAS does not provide for a 64-bit integer.  In
	    addition, SGI's older cblas.h can cause spurious warnings when
	    using the C-BLAS interface.

	1) most preferred:  use the optimized vendor-supplied library (such as
	    the Sun Performance Library, or IBM's ESSL).  This is often the
	    fastest, but might not be portable and might not always be
	    available.  When compiling a MATLAB mexFunction it might be
	    difficult get the mex compiler script to recognize the vendor-
	    supplied BLAS.

	2) When compiling the UMFPACK mexFunction to use UMFPACK in MATLAB, use
	    the BLAS provided by The Mathworks, Inc.  This assumes you are using
	    MATLAB V6 or higher, since the BLAS are not incorporated in V5 or
	    earlier versions.  On my Sun workstation, the MATLAB BLAS gave
	    slightly worse performance than the Sun Perf. BLAS.  The advantage
	    of using the MATLAB BLAS is that it's available on any computer that
	    has MATLAB V6 or higher.  I have not tried using MATLAB BLAS outside
	    of a mexFunction in a stand-alone C code, but MATLAB (V6) allows for
	    this.  This is well worth trying if you have MATLAB and don't want
	    to bother installing the ATLAS BLAS (option 3a, below).  The only
	    glitch to this is that MATLAB does not provide a portable interface
	    to the BLAS (an underscore is required for some but not all
	    architectures).  These variations are taken into account in the
	    mexopts.sh file provided with UMFPACK.  For Windows and MATLAB 6.0
	    or 6.1, you also need to copy the libmwlapack.dll file; see the
	    User Guide for details.

	    In the current distribution, the only BLAS that the UMFPACK 
	    mexFunction will use is the internal MATLAB BLAS.  It's possible to
	    use other BLAS, but handling the porting of using the mex compiler
	    with different BLAS libraries is not trivial.

	3) Use a portable high-performance BLAS library:

	    (a) The ATLAS BLAS, available at http://www.netlib.org/atlas,
		by R. Clint Whaley, Antoine Petitet, and Jack Dongarra.
		This has a standard C interface, and thus the interface to it is
		fully portable.  Its performance rivals, and sometimes exceeds,
		the vendor-supplied BLAS on many computers.

	    (b) The Fortran RISC BLAS by Michel Dayde', Iain Duff, Antoine
		Petitet, and Abderrahim Qrichi Aniba, available via anonymous
		ftp to ftp.enseeiht.fr in the pub/numerique/BLAS/RISC directory,
		See M. J. Dayde' and I. S. Duff, "The RISC BLAS:  A blocked
		implementation of level 3 BLAS for RISC processors, ACM Trans.
		Math. Software, vol. 25, no. 3., Sept. 1999.  This will give
		you good performance, but with the same C-to-Fortran portability
		problems as option (1).

	4) Use UMFPACK's built-in vanilla C code by setting -DNBLAS at compile
	    time.  The key advantage is portability, which is guaranteed if you
	    have an ANSI C compliant compiler.  You also don't need to download
	    any other package - UMFPACK is stand-alone.  No Fortran is used
	    anywhere in UMFPACK.  UMFPACK will be much slower than when using
	    options (1) through (3), however.

	5) least preferred:  use the standard Fortran implementation of the
	    BLAS, also available at Netlib (http://www.netlib.org/blas).  This
	    will be no faster than option (4), and not portable because of
	    C-to-Fortran calling conventions.  Don't bother trying option (5).

    The mechanics of how C calls the BLAS on various computers are as follows:

	* C-BLAS (from the ATLAS library, for example):
	    The same interface is used on all computers.

	* Defaults for calling the Fortran BLAS:
	    add underscore, pass scalars by reference, use string arguments.

	* The Fortran BLAS on Sun Solaris (when compiling the MATLAB mexFunction
	    or when using the Fortran RISC BLAS), SGI IRIX, Linux, and Compaq
	    Alpha: use defaults.

	* Sun Solaris (when using the C-callable Sun Performance library):
	    no underscore, pass scalars by value, use character arguments.

	* The Fortran BLAS (ESSL Library) on the IBM RS 6000, and HP Unix: 
	    no underscore, pass scalars by reference, use string arguments.

	* The Fortran BLAS on Windows:
	    no underscore, pass scalars by reference, use string arguments.
	    If you compile the umfpack mexFunction using umfpack_make, and are
	    using the lcc compiler bundled with MATLAB, then you must first
	    copy the umfpack\lcc_lib\libmwlapack.lib file into the
	    <matlab>\extern\lib\win32\lcc\ directory, where <matlab> is the
	    directory in which MATLAB is installed.  Next, type mex -setup
	    at the MATLAB prompt, and ask MATLAB to select the lcc compiler.
	    MATLAB has built-in BLAS, but it cannot be accessed by a program
	    compiled by lcc without first copying this file.
*/



/* -------------------------------------------------------------------------- */
#ifdef USE_C_BLAS	/* { */
/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
/* use the C-BLAS (any computer) */
/* -------------------------------------------------------------------------- */

/*
    C-BLAS is the default interface, with the following exceptions.  Solaris
    uses the Sun Performance BLAS for umfpack.a (the C-callable library).  SGI
    IRIX uses the SCSL BLAS for umfpack.a.  All architectures use MATLAB's
    internal BLAS for the mexFunction on any architecture.  These options are
    set in the Make.* files.  The Make.generic file uses no BLAS at all.

    If you use the ATLAS C-BLAS, then be sure to set the -I flag to
    -I/path/ATLAS/include, where /path/ATLAS is the ATLAS installation
    directory.  See Make.solaris for an example.  You do not need to do this
    for the SGI, which has a /usr/include/cblas.h.
*/

#include "cblas.h"

#ifdef COMPLEX
#define BLAS_GEMM_ROUTINE cblas_zgemm
#define BLAS_GEMV_ROUTINE cblas_zgemv
#define BLAS_GER_ROUTINE  cblas_zgeru
#define BLAS_DECLARE_SCALAR(x) double x [2]
#define BLAS_ASSIGN(x,xr,xi) { x [0] = xr ; x [1] = xi ; }
#else
#define BLAS_GEMM_ROUTINE cblas_dgemm
#define BLAS_GEMV_ROUTINE cblas_dgemv
#define BLAS_GER_ROUTINE  cblas_dger
#define BLAS_DECLARE_SCALAR(x) double x
#define BLAS_ASSIGN(x,xr,xi) { x = xr ; }
#endif

#define BLAS_NO_TRANSPOSE CblasNoTrans
#define BLAS_TRANSPOSE CblasTrans
#define BLAS_COLUMN_MAJOR_ORDER CblasColMajor,
#define BLAS_SCALAR(x) x
#define BLAS_INT_SCALAR(n) n
#define BLAS_ARRAY(a) a



/* -------------------------------------------------------------------------- */
#else	/* } USE_C_BLAS { */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* use Fortran (or other architecture-specific) BLAS */
/* -------------------------------------------------------------------------- */

/* No such argument when not using the C-BLAS */
#define BLAS_COLUMN_MAJOR_ORDER

/* Determine which architecture we're on and set options accordingly. */
/* The default, if nothing is defined is to add an underscore, */
/* pass scalars by reference, and use string arguments. */

/* ---------------------------------- */
/* Sun Performance BLAS */
/* ---------------------------------- */

#ifdef USE_SUNPERF_BLAS
#ifdef _SUNPERF_H
/* <sunperf.h> has been included somehow anyway, outside of umf_config.h */
#error "sunperf.h must NOT be #include'd.  See umf_config.h for details."
#endif
#define BLAS_BY_VALUE
#define BLAS_NO_UNDERSCORE
#define BLAS_CHAR_ARG
#endif	/* USE_SUNPERF_BLAS */

/* ---------------------------------- */
/* SGI SCSL BLAS */
/* ---------------------------------- */

#ifdef USE_SCSL_BLAS
#if defined (LP64)
#include <scsl_blas_i8.h>
#else
#include <scsl_blas.h>
#endif
#define BLAS_BY_VALUE
#define BLAS_NO_UNDERSCORE
#endif	/* USE_SCSL_BLAS */

/* ---------------------------------- */
/* IBM AIX, Windows, and HP Fortran BLAS */
/* ---------------------------------- */

#if defined (UMF_AIX) || defined (UMF_WINDOWS) || defined (UMF_HP)
#define BLAS_NO_UNDERSCORE
#endif


/* -------------------------------------------------------------------------- */
/* BLAS names */
/* -------------------------------------------------------------------------- */

#if defined (LP64) && defined (USE_SUNPERF_BLAS) && defined (LONG_INTEGER)

/* 64-bit sunperf BLAS, for Sun Solaris only */
#ifdef COMPLEX
#define BLAS_GEMM_ROUTINE zgemm_64
#define BLAS_GEMV_ROUTINE zgemv_64
#define BLAS_GER_ROUTINE  zgeru_64
#else
#define BLAS_GEMM_ROUTINE dgemm_64
#define BLAS_GEMV_ROUTINE dgemv_64
#define BLAS_GER_ROUTINE  dger_64
#endif	/* COMPLEX */

#else

#ifdef COMPLEX

/* naming convention (use underscore, or not) */
#ifdef BLAS_NO_UNDERSCORE
#define BLAS_GEMM_ROUTINE zgemm
#define BLAS_GEMV_ROUTINE zgemv
#define BLAS_GER_ROUTINE  zgeru
#else
/* default:  add underscore */
#define BLAS_GEMM_ROUTINE zgemm_
#define BLAS_GEMV_ROUTINE zgemv_
#define BLAS_GER_ROUTINE  zgeru_
#endif

#else

/* naming convention (use underscore, or not) */
#ifdef BLAS_NO_UNDERSCORE
#define BLAS_GEMM_ROUTINE dgemm
#define BLAS_GEMV_ROUTINE dgemv
#define BLAS_GER_ROUTINE  dger
#else
/* default:  add underscore */
#define BLAS_GEMM_ROUTINE dgemm_
#define BLAS_GEMV_ROUTINE dgemv_
#define BLAS_GER_ROUTINE  dger_
#endif

#endif	/* COMPLEX */

#endif /* LP64 && USE_SUNPERF_BLAS */


/* -------------------------------------------------------------------------- */
/* BLAS real or complex floating-point scalars */
/* -------------------------------------------------------------------------- */

#ifdef COMPLEX

/*
    The SunPerf BLAS expects to see a doublecomplex scalar, but it
    also will accept an array of size 2.  See the manual, normally at
    file:///opt/SUNWspro/WS6U1/lib/locale/C/html/manuals/perflib/user_guide
    /plug_using_perflib.html .  This manual is inconsistent with the man pages
    for zgemm, zgemv, and zgeru and also inconsistent with the <sunperf.h>
    include file.  Use this instead, for SunPerf (only works if you do NOT
    include sunperf.h).  Fortunately, this file (umf_config.h) is not included
    in any user code that calls UMFPACK.  Thus, the caller may include
    sunperf.h in his or her own code, and that is safely ignored here.
    SGI's SCSL BLAS has yet a different kind of struct, but we can use a
    double array of size 2 instead (since SCSL_VOID_ARGS is defined).
    Most BLAS expect complex scalars as pointers to double arrays of size 2.
*/

#define BLAS_DECLARE_SCALAR(x) double x [2]
#define BLAS_ASSIGN(x,xr,xi) { x [0] = xr ; x [1] = xi ; }
#define BLAS_SCALAR(x) x

#else

#define BLAS_DECLARE_SCALAR(x) double x
#define BLAS_ASSIGN(x,xr,xi) { x = xr ; }
#ifdef BLAS_BY_VALUE
#define BLAS_SCALAR(x) x
#else
#define BLAS_SCALAR(x) &(x)
#endif

#endif /* COMPLEX */


/* -------------------------------------------------------------------------- */
/* BLAS integer scalars */
/* -------------------------------------------------------------------------- */

/*
    Fortran requires integers to be passed by reference.
    The SCSL BLAS requires long long arguments in LP64 mode.
*/

#if defined (USE_SCSL_BLAS) && defined (LP64)
#define BLAS_INT_SCALAR(n) ((long long) n)
#else
#ifdef BLAS_BY_VALUE
#define BLAS_INT_SCALAR(n) n
#else
#define BLAS_INT_SCALAR(n) &(n)
#endif
#endif


/* -------------------------------------------------------------------------- */
/* BLAS strings */
/* -------------------------------------------------------------------------- */

/*
    The Sun Performance BLA wants a character instead of a string.
*/

#ifdef BLAS_CHAR_ARG
#define BLAS_NO_TRANSPOSE 'N'
#define BLAS_TRANSPOSE 'T'
#else
#define BLAS_NO_TRANSPOSE "N"
#define BLAS_TRANSPOSE "T"
#endif


/* -------------------------------------------------------------------------- */
/* BLAS arrays */
/* -------------------------------------------------------------------------- */

/*
    The complex SunPerf BLAS expects to see a doublecomplex array of size s.
    This is broken (see above, regarding complex scalars in sunperf.h).
    For SunPerf BLAS, just pass a pointer to the array, and ignore sunperf.h.
    With sunperf.h, you would need:

	#define BLAS_ARRAY(a) ((doublecomplex *)(a))

    SGI's SCSL BLAS has yet a different kind of struct, but we can use a
    double array of size 2 instead (since SCSL_VOID_ARGS is defined).

    The real versions all use just a (double *) pointer.

    In all cases, no typecast is required.  This will break if <sunperf.h> is
    included.

    If you have read this far, I hope you see now why (void *) a much better
    choice for complex BLAS prototypes, and why double x [2] is better than
    an architecture dependent struct { double real ; double imag ; }
    type definition.

*/

#define BLAS_ARRAY(a) (a)


/* -------------------------------------------------------------------------- */
#endif /* USE_C_BLAS } */
/* -------------------------------------------------------------------------- */





/* -------------------------------------------------------------------------- */
/* BLAS macros, for all interfaces */
/* -------------------------------------------------------------------------- */

/*
   All architecture dependent issues have now been taken into consideration,
   and folded into the macros BLAS_DECLARE_SCALAR, BLAS_ASSIGN, BLAS_*_ROUTINE,
   BLAS_COLUMN_MAJOR_ORDER, BLAS_NO_TRANSPOSE, BLAS_TRANSPOSE, BLAS_SCALAR,
   BLAS_INT_SCALAR, BLAS_ARRAY, and Int.

   You will note that there is not a *** single *** name, declaration, or
   argument to the BLAS which is not somehow different in one or more versions
   of the BLAS!
*/




/* C = C - A*B, where A is m-by-k, B is k-by-n, and leading dimension is d */
#define BLAS_GEMM(m,n,k,A,B,C,d) \
{ \
    BLAS_DECLARE_SCALAR (alpha) ; \
    BLAS_DECLARE_SCALAR (beta) ; \
    BLAS_ASSIGN (alpha, -1.0, 0.0) ; \
    BLAS_ASSIGN (beta, 1.0, 0.0) ; \
    (void) BLAS_GEMM_ROUTINE (BLAS_COLUMN_MAJOR_ORDER \
	BLAS_NO_TRANSPOSE, BLAS_NO_TRANSPOSE, \
	BLAS_INT_SCALAR (m), BLAS_INT_SCALAR (n), BLAS_INT_SCALAR (k), \
	BLAS_SCALAR (alpha), \
	BLAS_ARRAY (A), BLAS_INT_SCALAR (d), \
	BLAS_ARRAY (B), BLAS_INT_SCALAR (d), BLAS_SCALAR (beta), \
	BLAS_ARRAY (C), BLAS_INT_SCALAR (d)) ; \
}


/* A = A - x*y', where A is m-by-n with leading dimension d */
#define BLAS_GER(m,n,x,y,A,d) \
{ \
    Int incx = 1 ; \
    BLAS_DECLARE_SCALAR (alpha) ; \
    BLAS_ASSIGN (alpha, -1.0, 0.0) ; \
    (void) BLAS_GER_ROUTINE (BLAS_COLUMN_MAJOR_ORDER \
	BLAS_INT_SCALAR (m), BLAS_INT_SCALAR (n), \
	BLAS_SCALAR (alpha), \
	BLAS_ARRAY (x), BLAS_INT_SCALAR (incx), \
	BLAS_ARRAY (y), BLAS_INT_SCALAR (d), \
	BLAS_ARRAY (A), BLAS_INT_SCALAR (d)) ; \
}


/* y = y - A'*x, where A is m-by-n with leading dimension d, */
/* and x and y are row vectors with stride d */
#define BLAS_GEMV_ROW(m,n,A,x,y,d) \
{ \
    BLAS_DECLARE_SCALAR (alpha) ; \
    BLAS_DECLARE_SCALAR (beta) ; \
    BLAS_ASSIGN (alpha, -1.0, 0.0) ; \
    BLAS_ASSIGN (beta, 1.0, 0.0) ; \
    (void) BLAS_GEMV_ROUTINE (BLAS_COLUMN_MAJOR_ORDER \
	BLAS_TRANSPOSE, \
	BLAS_INT_SCALAR (m), BLAS_INT_SCALAR (n), \
	BLAS_SCALAR (alpha), \
	BLAS_ARRAY (A), BLAS_INT_SCALAR (d), \
	BLAS_ARRAY (x), BLAS_INT_SCALAR (d), BLAS_SCALAR (beta), \
	BLAS_ARRAY (y), BLAS_INT_SCALAR (d)) ; \
}


/* y = y - A*x, where A is m-by-n with leading dimension d, */
/* and x and y are column vectors with stride 1 */
#define BLAS_GEMV_COL(m,n,A,x,y,d) \
{ \
    Int incx = 1 ; \
    Int incy = 1 ; \
    BLAS_DECLARE_SCALAR (alpha) ; \
    BLAS_DECLARE_SCALAR (beta) ; \
    BLAS_ASSIGN (alpha, -1.0, 0.0) ; \
    BLAS_ASSIGN (beta, 1.0, 0.0) ; \
    (void) BLAS_GEMV_ROUTINE (BLAS_COLUMN_MAJOR_ORDER \
	BLAS_NO_TRANSPOSE, \
	BLAS_INT_SCALAR (m), BLAS_INT_SCALAR (n), \
	BLAS_SCALAR (alpha), \
	BLAS_ARRAY (A), BLAS_INT_SCALAR (d), \
	BLAS_ARRAY (x), BLAS_INT_SCALAR (incx), BLAS_SCALAR (beta), \
	BLAS_ARRAY (y), BLAS_INT_SCALAR (incy)) ; \
}


#endif	/* !defined (USE_NO_BLAS) } */



