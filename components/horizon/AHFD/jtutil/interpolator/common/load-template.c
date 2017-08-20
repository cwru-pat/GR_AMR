/* load-template.c -- template for load functions */
/* $Header$ */

/*
 * Each of the functions defined in this file loades a molecule-sized
 * piece of an input array into a data struct.  There is one function for
 * each combination of molecule size and input array datatype.
 */

/******************************************************************************/

/*
 * load-data routines for real datatypes
 */

#undef  DATA
#define DATA	DATA_REAL

void LOAD_FUNCTION_NAME(r)(const CCTK_REAL *ptr,
			   INT_STRIDE_IJK,
			  struct DATA_STRUCT *data)
{
#include LOAD_DATA_FILE_NAME
}

#ifdef HAVE_CCTK_REAL4
  void LOAD_FUNCTION_NAME(r4)(const CCTK_REAL4 *ptr,
			      INT_STRIDE_IJK,
			      struct DATA_STRUCT *data)
  {
  #include LOAD_DATA_FILE_NAME
  }
#endif

#ifdef HAVE_CCTK_REAL8
  void LOAD_FUNCTION_NAME(r8)(const CCTK_REAL8 *ptr,
			      INT_STRIDE_IJK,
			      struct DATA_STRUCT *data)
  {
  #include LOAD_DATA_FILE_NAME
  }
#endif

#ifdef HAVE_CCTK_REAL16
  void LOAD_FUNCTION_NAME(r16)(const CCTK_REAL16 *ptr,
				INT_STRIDE_IJK,
				struct DATA_STRUCT *data)
  {
  #include LOAD_DATA_FILE_NAME
  }
#endif

/******************************************************************************/

/*
 * load-data routines for complex datatypes
 */

#undef  DATA
#define DATA	DATA_COMPLEX

void LOAD_FUNCTION_NAME(c)(const CCTK_REAL (*ptr)[COMPLEX_N_PARTS],
			   INT_STRIDE_IJK, int part,
			   struct DATA_STRUCT *data)
{
#include LOAD_DATA_FILE_NAME
}

#ifdef HAVE_CCTK_COMPLEX8
  void LOAD_FUNCTION_NAME(c8)(const CCTK_REAL4 (*ptr)[COMPLEX_N_PARTS],
			      INT_STRIDE_IJK, int part,
			      struct DATA_STRUCT *data)
  {
  #include LOAD_DATA_FILE_NAME
  }
#endif

#ifdef HAVE_CCTK_COMPLEX16
  void LOAD_FUNCTION_NAME(c16)(const CCTK_REAL8 (*ptr)[COMPLEX_N_PARTS],
				INT_STRIDE_IJK, int part,
				struct DATA_STRUCT *data)
  {
  #include LOAD_DATA_FILE_NAME
  }
#endif

#ifdef HAVE_CCTK_COMPLEX32
  void LOAD_FUNCTION_NAME(c32)(const CCTK_REAL16 (*ptr)[COMPLEX_N_PARTS],
				INT_STRIDE_IJK, int part,
				struct DATA_STRUCT *data)
  {
  #include LOAD_DATA_FILE_NAME
  }
#endif
