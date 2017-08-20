/* lapack.h -- C/C++ prototypes for (some) BLAS+LAPACK+wrapper routines */
/* $Header$ */

/*
 * prerequisites:
 *	"cctk.h"
 *	"config.h"	// for "integer" = Fortran integer
 */

#ifdef __cplusplus
extern "C"
       {
#endif

/*
 * ***** BLAS *****
 */
integer CCTK_FCALL
  CCTK_FNAME(isamax)(const integer* N, const float SX[], const integer* incx);
integer CCTK_FCALL
  CCTK_FNAME(idamax)(const integer* N, const double DX[], const integer* incx);

/*
 * ***** LAPACK *****
 */
void CCTK_FCALL
  CCTK_FNAME(sgesv)(const integer* N, const integer* NRHS,
		    float A[], const int* LDA,
		    integer IPIV[],
		    float B[], const integer* LDB, integer* info);
void CCTK_FCALL
  CCTK_FNAME(dgesv)(const integer* N, const integer* NRHS,
		    double A[], const int* LDA,
		    integer IPIV[],
		    double B[], const integer* LDB, integer* info);

/*
 * ***** wrappers (for passing character-string args) *****
 */
/* norm_int = 0 for infinity-norm, 1 for 1-norm */
void CCTK_FCALL
  CCTK_FNAME(sgecon_wrapper)(const integer* norm_int,
			     const integer* N,
			     float A[], const integer* LDA,
			     const float* anorm, float* rcond,
			     float WORK[], integer IWORK[],
			     integer* info);
void CCTK_FCALL
  CCTK_FNAME(dgecon_wrapper)(const integer* norm_int,
			     const integer* N,
			     double A[], const integer* LDA,
			     const double* anorm, double* rcond,
			     double WORK[], integer IWORK[],
			     integer* info);

#ifdef __cplusplus
       }	/* extern "C" */
#endif
