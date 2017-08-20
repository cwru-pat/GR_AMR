 /*@@
   @file      util.c
   @date      23 Oct 2001
   @author    Jonathan Thornburg <jthorn@aei.mpg.de>
   @desc
	Utility routines for generalized interpolation.

	This file contains various utility routines for the interpolator.
   @enddesc
   @version   $Id$
 @@*/

#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../../AHFD_macros.h"

#include "InterpLocalUniform.h"

/* the rcs ID and its dummy function to use it */
/* static const char *rcsid = "$Header$"; */
/* #ifndef AEILOCALINTERP_STANDALONE_TEST */
/*   CCTK_FILEVERSION(AEIThorns_AEILocalInterp_src_util_c) */
/* #endif */

/******************************************************************************/

/*@@
  @routine  AEILocalInterp_decode_N_parts
  @date     22 Jan 2002
  @author   Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc     This function decodes a CCTK_VARIABLE_* variable type code
	    (cctk_Constants.h) to determine whether the type is real
	    or complex.
  @enddesc

  @var      type_code
  @vdesc    The type code to be decoded
  @vtype    int
  @endvar

  @returntype	int
  @returndesc	This function returns
		1	if the data type represents a real number of some sort
			(includes integers)
		2	if the data type represents a complex number
		0	if the data type doesn't represent a number,
			eg strings and pointers
		-1	if the data type is invalid
  @endreturndesc
  @@*/
int AEILocalInterp_decode_N_parts(int type_code)
{
switch	(type_code)
	{
case CCTK_VARIABLE_VOID:	return 0;
case CCTK_VARIABLE_BYTE:	return 1;
case CCTK_VARIABLE_INT:		return 1;
case CCTK_VARIABLE_INT2:	return 1;
case CCTK_VARIABLE_INT4:	return 1;
case CCTK_VARIABLE_INT8:	return 1;
case CCTK_VARIABLE_REAL:	return 1;
case CCTK_VARIABLE_REAL4:	return 1;
case CCTK_VARIABLE_REAL8:	return 1;
case CCTK_VARIABLE_REAL16:	return 1;
case CCTK_VARIABLE_COMPLEX:	return 2;
case CCTK_VARIABLE_COMPLEX8:	return 2;
case CCTK_VARIABLE_COMPLEX16:	return 2;
case CCTK_VARIABLE_COMPLEX32:	return 2;
case CCTK_VARIABLE_STRING:	return 0;
case CCTK_VARIABLE_POINTER:	return 0;
case CCTK_VARIABLE_FPOINTER:	return 0;
default:			return -1;
	}
}

/******************************************************************************/

/*@@
  @routine  AEILocalInterp_get_int_param
  @date     12 Jan 2007
  @author   Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc     This function gets the (scalar) value of a CCTK_INT or Boolean
            Cactus parameter.
  @enddesc

  @var      thorn_or_implementation_name
  @vdesc    The name of the thorn (for private parameters) or implementation
            (for restricted parameters), eg "AEILocalInterp".
  @vtype    const char*
  @endvar

  @var      parameter_name
  @vdesc    The name of the parameter, e.g. "log_interp_coords"
  @vtype    const char*
  @endvar

  @returntype	int
  @returndesc	This function returns the value of the parameter.
                If an error occurs, this function calls
                  CCTK_VWarn(CCTK_WARN_ABORT, ...)
                (and does not return to the caller).
  @endreturndesc
  @@*/
int AEILocalInterp_get_int_param(const char* const thorn_or_implementation_name,
				 const char* const parameter_name)
{
  printf("%s","!!!!!! AEILocalInterp_get_int_param\n");
/* const CCTK_INT* const value_ptr */
/* 	= (const CCTK_INT*) CCTK_ParameterGet(parameter_name, */
/* 					      thorn_or_implementation_name, */
/*                                               NULL); */

/* if (value_ptr == NULL) */
/*    then CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING, */
/* "***** AEILocalInterp_decode_N_parts():\n" */
/* "        can't get value of parameter %s::%s!\n" */
/* 		   , */
/* 		   thorn_or_implementation_name, parameter_name); /\*NOTREACHED*\/ */

/* return *value_ptr; */
  return 0;
}
