/* ilucg.h -- C/C++ prototypes for [sd]ilucg() routines */
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
 * ***** ILUCG *****
 */
void CCTK_FCALL
  CCTK_FNAME(silucg)(const integer* N,
		     const integer IA[], const integer JA[], const float A[],
		     const float B[], float X[],
		     integer ITEMP[], float RTEMP[],
		     const float* EPS, const integer* ITER,
		     integer* ISTATUS);
void CCTK_FCALL
  CCTK_FNAME(dilucg)(const integer* N,
		     const integer IA[], const integer JA[], const double A[],
		     const double B[], double X[],
		     integer ITEMP[], double RTEMP[],
		     const double* EPS, const integer* ITER,
		     integer* ISTATUS);

#ifdef __cplusplus
       }	/* extern "C" */
#endif
