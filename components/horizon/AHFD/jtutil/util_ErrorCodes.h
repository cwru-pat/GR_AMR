/*@@
  @file         util_ErrorCodes.h
  @date         Wed Nov  7 16:16:08 CET 2001
  @author       Jonathan Thornburg <jthorn@aei.mpg.de>
  @header       util_ErrorCodes.h
  @version      $Header$
   @desc 
        This header defines the error codes returned by Util_* functions.
   @enddesc 
  @@*/

#ifndef _UTIL_ERRORCODES_H_
#define _UTIL_ERRORCODES_H_ 1

/******************************************************************************/

/*
 * These are generic error codes, used by multiple Util_*() functions.
 * They are between -1 and -99 inclusive.  (All Cactus error codes are
 * -ve integers.)
 */

/*@@
  @defines      UTIL_ERROR_NO_MEMORY
  @desc         error return code: unable to allocate memory
  @@*/
#define UTIL_ERROR_NO_MEMORY    (-1)

/*@@
  @defines      UTIL_ERROR_BAD_HANDLE
  @desc         error return code: handle is invalid
  @@*/
#define UTIL_ERROR_BAD_HANDLE   (-2)

/*@@
  @defines      UTIL_ERROR_BAD_INPUT
  @desc         error return code: input is invalid
  @@*/
#define UTIL_ERROR_BAD_INPUT    (-3)

/******************************************************************************/

#endif  /* _UTIL_ERRORCODES_H_ */
