#ifndef COSMO_AHFD_MACROS
#define COSMO_AHFD_MACROS

#define then	/* empty */


#define AHFD AHFinderDirect

#define CCTK_FCALL

#define CCTK_FNAME(fun) fun##_


#define HAVE_CAPABILITY_HDF5 HAVE_HDF5

//#define __cplusplus
typedef int    CCTK_INT;
typedef double CCTK_REAL;
typedef int integer;

typedef CCTK_REAL fp;

typedef void *CCTK_POINTER;
typedef const void *CCTK_POINTER_TO_CONST;
typedef void (*CCTK_FPOINTER)(void);
typedef unsigned char CCTK_BYTE;

#define HAVE_CCTK_POINTER 1
#define HAVE_CCTK_POINTER_TO_CONST 1
#define HAVE_CCTK_FPOINTER 1

/* Character types */
typedef char CCTK_CHAR;
typedef const char * CCTK_STRING;
#define HAVE_CCTK_CHAR 1
#define HAVE_CCTK_STRING 1

#define CCODE

#define FATAL_ERROR	(-1)
#define FP_IS_DOUBLE

#define FP_SCANF_FORMAT "%lf"
#define FP_EPSILON	DBL_EPSILON	/* from <float.h> */



#define FINITE_DIFF_ORDER	4

/* store as (Fortran) dense matrix, solve with LAPACK */
#undef  HAVE_DENSE_JACOBIAN__LAPACK

/* store as row-oriented sparse matrix, solve with ILUCG */
#define HAVE_ROW_SPARSE_JACOBIAN__ILUCG

/* store as row-oriented sparse matrix, solve with UMFPACK */
#undef HAVE_ROW_SPARSE_JACOBIAN__UMFPACK

#define FORTRAN_INTEGER_IS_INT
#undef  FORTRAN_INTEGER_IS_LONG

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*
 * The following #ifdefs and #defines shouldn't need changing; they
 * decode the previous set of Jacobian storage format and linear solver
 * #defines to determine which storage formats we are using.
 */
#ifdef HAVE_DENSE_JACOBIAN__LAPACK
#define HAVE_DENSE_JACOBIAN
#endif

#ifdef HAVE_ROW_SPARSE_JACOBIAN__ILUCG
#define HAVE_ROW_SPARSE_JACOBIAN
#endif

#ifdef HAVE_ROW_SPARSE_JACOBIAN__UMFPACK
#define HAVE_ROW_SPARSE_JACOBIAN
#endif

/******************************************************************************/

/*
 * misc stuff
 */


#define STRING_EQUAL(s_,t_)	(strcmp(s_,t_) == 0)

/* C library uses abs() for integer absolute value; I prefer iabs() */
#define iabs(x_)	abs(x_)

#define MAXDIM 3

#define RETURN_SUCCESS 1
/******************************************************************************/

/*
 * misc sentinal values
 */

/* n.b. this value is usually 10...0 binary for 2's-complement arithmetic */
#define INT_SENTINAL_VALUE	(~ INT_MAX)	/* from <limits.h> */

#define FLOAT_SENTINAL_VALUE	(- FLT_MAX)	/* from <float.h> */
#define DOUBLE_SENTINAL_VALUE	(- DBL_MAX)	/* from <float.h> */

#define error_exit(a, ...)                       \
  do{                                              \
    char s[200];                                   \
    snprintf(s, 200, __VA_ARGS__);                 \
    TBOX_ERROR(s);                                 \
   }while(0)


#define CCTK_VError(a, b, c, ...)                \
  error_exit(c, __VA_ARGS__)

#define CCTK_VWarn(a, b, c, d, ...)                \
  do{                                              \
    char s[200];                                   \
    snprintf(s, 200, __VA_ARGS__);                 \
    TBOX_WARNING(s);                               \
   }while(0)

#define CCTK_CVWarn(a, b, c, d, ...)                \
  do{                                              \
    char s[200];                                   \
    snprintf(s, 200, __VA_ARGS__);                 \
    printf("%s",s);                                \
   }while(0)


#define CCTK_WARN(a, b) \
  tbox::pout<<"Waning!!!\n"<<b<<"\n"

#define CCTK_VInfo(a, ...) \
  do{                                           \
    char s[200];                                \
    snprintf(s, 200, __VA_ARGS__);              \
    tbox::pout<<s<<"\n";                        \
  }while(0)
//printf( __VA_ARGS__), std::cout<<"\n"



  
/* #define assert(a)                               \ */
/*   TBOX_ASSERT(a) */


#define PASS_GROUPSIZE(group, dir)  CCTKGROUPNUM_##group >= 0 ?         \
                                    CCTK_ArrayGroupSizeI(GH, dir, CCTKGROUPNUM_##group) : &_cctk_zero

#define PASS_GROUPLEN(thorn, group) CCTKGROUPNUM_##group >= 0 ? \
                                    CCTKi_GroupLengthAsPointer(#thorn "::" #group) : &_cctk_zero

#define PASS_REFERENCE(var, level)  CCTKARGNUM_##var >= 0 ? \
                                    GH->data[CCTKARGNUM_##var][level] : 0

#define CCTK_ARGUMENTS AHFINDERDIRECT_CARGUMENTS
#define _CCTK_ARGUMENTS _CCTK_CARGUMENTS
#define DECLARE_CCTK_ARGUMENTS DECLARE_AHFINDERDIRECT_CARGUMENTS

  
/* #define Util_TableQueryValueInfo(a, b, c, d) \ */
/*   db->keyExists(d) */


/* #define Util_TableGetString(a, b, c, d) \ */
/*   db->keyExist(d); if(db->keyExist(d)) c = db->getString(d).c_str() */

/* #define Util_TableSetString(a, b, c) \ */
/*    - (int)db->keyExist(c); if(!db->keyExist(c)) db->putString(c, toStr(b));  */
  
/*@@
  @routine INTERPOLATE (macro)
  @date    18 Oct 2001
  @author  code by ???, these comments by Jonathan Thornburg
  @desc
       This macro does the interpolation of in_array[] to compute
       a single value  out_array[n]  (actually  out_array[n]subpart ;
       see the comments below for details).  The data to be interpolated
       must be real numbers of some type, i.e. if the arrays are
       complex this macro must be called separately for the real
       and imaginary parts.
  @enddesc

  @var    cctk_type
  @vdesc  C type of input and output array elements (might be complex)
  @endvar

  @var    cctk_subtype
  @vdesc  C type of actual numbers being interpolated (must be real)
  @endvar

  @var    subpart
  @vdesc  string to be suffixed to input/output array element to get
          to get real number, i.e. empty string if cctk_type is real,
          .Re or .Im as appropriate if cctk_type is complex
  @endvar

  @var    in_array
  @vdesc  A pointer to array to be interpolated (strictly speaking, to
          the array's [0][0]...[0] element); this is typically passed
          as a  void *  pointer so we typecast it as necessary.
  @endvar

  @var    out_array
  @vdesc  A 1-dimensional array where the interpolation result should be
          stored; this is typically passed as a  void *  pointer so we
          typecast it as necessary.
  @endvar

  @var    order
  @vdesc  The order of the interpolation (1=linear, 2=quadratic, 3=cubic, ...)
  @endvar

  @var    point
  @vdesc  [MAXDIM] array of integers giving the integer grid coordinates
          of the closest grid point to the interpolation point; the
          interpolation stencil/molecule is centered at this point.
  @endvar

  @var    dims
  @vdesc  [MAXDIM] array of integers giving the dimensions of  in_array .
  @endvar

  @var    n
  @vdesc  Position in  out_array  where we should store the interpolation
          result.
  @endvar

  @var    coeff
  @vdesc  [MAXDIM][MAX_ORDER+1] array of (floating-point) interpolation
          coefficients; detailed semantics are that coeff[axis][m] is the
          coefficient of y[m] when the 1-dimensional Lagrange interpolation
          polynomial passing through the  order+1  points
             {(0,y[0]), (1,y[1]), ..., (order,y[order])}
          is evaluated at the position x=offset[axis].
  @endvar
  @@*/
/*
 * The basic idea here is that conceptually we first interpolate the
 * (say) 3D gridfn in the x direction at each y and z grid point,
 * then interpolate that 2D plane of values in the y direction at
 * each z grid point, and finally interpolate that 1D line of values
 * in the z direction.  The implementation actually interleaves the
 * different directions' interpolations so that only 3 scalar temporaries
 * are needed.
 */
#define INTERPOLATE(cctk_type, in_array, out_array,                           \
                    order, point, dims, n, coeff)                             \
    {                                                                         \
      int ii, jj, kk;                                                         \
      const cctk_type *fi;                                                    \
      cctk_type interp_result, fj, fk;                                        \
                                                                              \
                                                                              \
      interp_result = 0;                                                      \
                                                                              \
      /* NOTE-MAXDIM: support >3D arrays by adding more loops */              \
      for (kk = 0; kk <= order; kk++)                                         \
      {                                                                       \
        fk = 0;                                                               \
        for (jj = 0; jj <= order; jj++)                                       \
        {                                                                     \
          /* NOTE-MAXDIM: for >3D arrays adapt the index calculation here */  \
          fi = (const cctk_type *) in_array +                                 \
               point[0] + dims[0]*(point[1]+jj + dims[1]*(point[2]+kk));      \
                                                                              \
          fj = 0;                                                             \
          for (ii = 0; ii <= order; ii++)                                     \
          {                                                                   \
            fj += fi[ii] * coeff[0][ii];                                      \
          }                                                                   \
          /* at this point we have just computed */                           \
          /* fj = in_array[*][jj][kk] interpolated to x=offset[0] */          \
                                                                              \
          fk += fj * coeff[1][jj];                                            \
        }                                                                     \
        /* at this point we have just computed */                             \
        /* fk = fj[*][kk] interpolated to y=offset[1] */                      \
        /*       = in_array[*][*][kk] interpolated to */                      \
        /*                            x=offset[0], y=offset[1] */             \
                                                                              \
        interp_result += fk * coeff[2][kk];                                   \
      }                                                                       \
      /* at this point we have just computed */                               \
      /* interp_result = fk[*] interpolated to z=offset[2] */                 \
      /*               = in_array[*][*][*] interpolated to */                 \
      /*                 x=offset[0], y=offset[1], z=offset[2] */             \
                                                                              \
      /* assign the result */                                                 \
      ((cctk_type *) out_array)[n] = interp_result;                           \
    }                                                      /* end of macro */

#define CCTK_VARIABLE_VOID             100
#define CCTK_VARIABLE_BYTE             110
#define CCTK_VARIABLE_INT              120
#define CCTK_VARIABLE_INT1             121
#define CCTK_VARIABLE_INT2             122
#define CCTK_VARIABLE_INT4             123
#define CCTK_VARIABLE_INT8             124
#define CCTK_VARIABLE_INT16            125
#define CCTK_VARIABLE_REAL             130
#define CCTK_VARIABLE_REAL4            131
#define CCTK_VARIABLE_REAL8            132
#define CCTK_VARIABLE_REAL16           133
#define CCTK_VARIABLE_COMPLEX          140
#define CCTK_VARIABLE_COMPLEX8         141
#define CCTK_VARIABLE_COMPLEX16        142
#define CCTK_VARIABLE_COMPLEX32        143
#define CCTK_VARIABLE_CHAR             150
#define CCTK_VARIABLE_STRING           151
#define CCTK_VARIABLE_POINTER          160
#define CCTK_VARIABLE_POINTER_TO_CONST 161
#define CCTK_VARIABLE_FPOINTER         162

/* DEPRECATED IN BETA 12 */
#define CCTK_VARIABLE_FN_POINTER CCTK_VARIABLE_FPOINTER

/* steerable status of parameters */
#define CCTK_STEERABLE_NEVER   200
#define CCTK_STEERABLE_ALWAYS  201
#define CCTK_STEERABLE_RECOVER 202

/* group distributions */
#define CCTK_DISTRIB_CONSTANT 301
#define CCTK_DISTRIB_DEFAULT  302

/* group types */
#define CCTK_SCALAR 401
#define CCTK_GF     402
#define CCTK_ARRAY  403

/* group scopes */
#define CCTK_PRIVATE   501
#define CCTK_PROTECTED 502
#define CCTK_PUBLIC    503

/* constants for CCTK_TraverseString() */
#define CCTK_VAR          601
#define CCTK_GROUP        602
#define CCTK_GROUP_OR_VAR 603

/* the grid is too small for the selected interpolation molecule */
#define CCTK_ERROR_INTERP_GRID_TOO_SMALL (-1000)
/* ... old code for backwards compatability */
#define CCTK_ERROR_INTERP_GRID_TOO_TINY  CCTK_ERROR_INTERP_GRID_TOO_SMALL

/*
 * the (multiprocessor) grid's ghostzone size is too small for the selected
 * interpolation molecule (or this processor's chunk of the grid is too small)
 */
#define CCTK_ERROR_INTERP_GHOST_SIZE_TOO_SMALL (-1001)

/*
 * an interpolation point is outside (or too close to an edge of)
 * the input grid
 */
#define CCTK_ERROR_INTERP_POINT_OUTSIDE (-1002)
/* ... old code for backwards compatability */
#define CCTK_ERROR_INTERP_POINT_X_RANGE CCTK_ERROR_INTERP_POINT_OUTSIDE

/* an interpolation point is in (or too close to) an excised region */
#define CCTK_ERROR_INTERP_POINT_EXCISED (-1003)

/* an interpolation coordinate (or some other intermediate value in the */
/* interpolation process) was an IEEE NaN or other non-finite number */
#define CCTK_ERROR_INTERP_COORD_NAN	(-1004)

/* the grid spacing was specified as zero along at least one axis */
#define CCTK_ERROR_INTERP_DELTA_X_ZERO	(-1005)

#ifndef S_ISDIR
#define S_ISDIR(mode)   (((mode) & S_IFMT) == S_IFDIR)
#endif

#define HAVE_MODE_T 1

/* Some systems don't have mode_t and only pass one argument to mkdir. */
#ifdef HAVE_MODE_T
#define MKDIR_WRAPPER(a,b) mkdir(a,b)
#else
#define MKDIR_WRAPPER(a,b) mkdir(a)
#endif

#define USE_HERMITE_GEOM_INTERP 1

#define DERIV(x)   x

#endif
