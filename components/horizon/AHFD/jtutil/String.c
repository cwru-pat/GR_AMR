 /*@@
   @file      String.c
   @date      Tue May  2 10:44:19 2000
   @author    Tom Goodale
   @desc
              Routines dealing with strings.
   @enddesc
   @version   $Id$
 @@*/

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include "util_String.h"

static const char *rcsid = "$Header$";


/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    Util_StrSep
   @date       Tue May  2 10:29:07 2000
   @author     Tom Goodale
   @desc
     The strsep() function returns the next token from the string stringp which is delimited by delim.  The token
     is terminated with a `\0' character and stringp is updated to point past the token.

     RETURN VALUE
     The strsep() function returns a pointer to the token, or NULL if delim is not found in stringp.

   @enddesc
   @calls
   @calledby
   @history

   @endhistory
   @var     stringp
   @vdesc   The string to search for a token in.
   @vtype   const char **stringp
   @vio     inout
   @vcomment

   @endvar
   @var     delim
   @vdesc   The delimiter
   @vtype   const char *delim
   @vio     in
   @vcomment

   @endvar

   @returntype const char *
   @returndesc
   a pointer to the token, or NULL if delim is not found in stringp.
   @endreturndesc

@@*/
const char *Util_StrSep(const char **stringp, const char *delim)
{
  int retlength = 0;
  static char *retval = NULL;
  char *temp;
  const char *start;
  const char *end;

  start = *stringp;

  end = strstr(start, delim);

  /* Is the delimiter part of the string */
  if(end)
  {
    if(retlength < (end-start)+1)
    {
      temp = realloc(retval, (end-start+1));

      if(temp)
      {
        retval = temp;
        retlength = end-start+1;
      }
      else
      {
        free(retval);
        retval = NULL;
        retlength = 0;
      }
    }

    if(retval)
    {
      strncpy(retval, start, (size_t)(end-start));
      retval[end-start] = '\0';

      *stringp = end+strlen(delim);
    }

  }
  else
  {
    free(retval);
    retval = NULL;
    retlength = 0;
  }

  return retval;
}

 /*@@
   @routine    Util_SplitString
   @date       Wed Jan 20 10:14:00 1999
   @author     Tom Goodale
   @desc
   Splits a string into two parts at the given seperator.
   Assigns memory for the two resulting strings, so this should be freed
   when no longer needed.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory
   @var     before
   @vdesc   String before seperator
   @vtype   char **
   @vio     out
   @vcomment

   @endvar
   @var     after
   @vdesc   String after seperator
   @vtype   char **
   @vio     out
   @vcomment

   @endvar
   @var     string
   @vdesc   String to seperate
   @vtype   const char *
   @vio     in
   @vcomment

   @endvar
   @var     sep
   @vdesc   String seperator
   @vtype   const char *
   @vio     in
   @vcomment

   @endvar

   @returntype int
   @returndesc
   0 - success
   1 - seperator not found
   2 - out of memory
   @endreturndesc

@@*/
int Util_SplitString(char **before, char **after, const char *string, const char *sep)
{
  int retval=0;
  char *position;

  /* Find location of the seperator */
  position = strstr(string, sep);

  if(position)
  {
    /*Allocate memory for return strings. */
    *before = (char *)malloc((size_t)((position-string+1)*sizeof(char)));
    *after  = (char *)malloc((size_t)((strlen(string)-(position-string)-strlen(sep)+1)*sizeof(char)));

    /* Check that the allocation succeeded. */
    if(!*before || !*after)
    {
      free(*before);
      *before = NULL;
      free(*after);
      *after = NULL;
      retval = 2;
    }
    else
    {
      retval = 3;
    }
  }
  else
  {
    *before = NULL;
    *after = NULL;
    retval = 1;
  }

  if(position && *before && *after)
  {
    /* Copy the data */
    strncpy(*before, string, (int)(position-string));
    (*before)[(int)(position-string)] = '\0';

    strncpy(*after, position+strlen(sep), strlen(string)-(int)(position-string)-strlen(sep));
    (*after)[strlen(string)-(position-string)-strlen(sep)] = '\0';

    retval = 0;
  }

  return retval;
}


/*@@
  @routine  Util_Strlcpy
  @date     1.Feb.2003
  @author   Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc     This function implements the  strlcpy()  function described in
               http://www.openbsd.org/papers/strlcpy-paper.ps

            The strlcpy(3) function copies up to  size-1  characters
            from the null-terminated string  src  to  dst , followed
            by a null character (so  dst  is always null-terminated).

            The strlcpy(3) API is intended to replace strncpy(3).  In
            comparison to strncpy(3), strlcpy(3) is safer and easier to
            use (it guarantees null termination of the destination buffer),
            and faster (it doesn't have to fill the entire buffer with
            null characters).
  @enddesc

  @var      dst
  @vdesc    A non-null pointer to the destination buffer.
  @vtype    char * 
  @endvar

  @var      src
  @vdesc    A non-null pointer to the source string.
  @vtype    const char *
  @endvar

  @var      dst_size
  @vdesc    The size of the destination buffer.
  @vtype    size_t
  @endvar

  @returntype   size_t
  @returndesc   This function returns strlen(src).
  @endreturndesc
@@*/
size_t Util_Strlcpy(char* dst, const char* src, size_t dst_size)
{
  const size_t src_size = strlen(src);
  if (src_size < dst_size)
  {
    strcpy(dst, src);
  }
  else
  {
    strncpy(dst, src, dst_size-1);
    dst[dst_size-1] = '\0';
  }
  return src_size;
}

/*@@
  @routine  Util_Strlcat
  @date     16.Feb.2003
  @author   Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc     This function implements the  strcat()  function described in
               http://www.openbsd.org/papers/strlcpy-paper.ps

            The strlcat(3) function appends the null-terminated string
            src to the end of dst.  It will append at most
                size - strlen(dst) - 1
            characters, and always null-terminates the result.
            (Hence this function never overflows the destination buffer.)

            The strlcat(3) is intended to replace strncat(3).  In
            comparison to strncat(3), strlcat(3) is safer and easier to
            use: it guarantees null termination of the destination buffer,
            and the size parameter is easy to specify without danger
            of off-by-one errors.
  @enddesc

  @var      dst
  @vdesc    A non-null pointer to the destination buffer.
  @vtype    char *
  @endvar

  @var      src
  @vdesc    A non-null pointer to the source string.
  @vtype    const char *
  @endvar

  @var      dst_size
  @vdesc    The size of the destination buffer.
  @vtype    size_t
  @endvar

  @returntype   size_t
  @returndesc   This function returns the length of the string it
                tries to create, i.e.  strlen(src) + strlen(dst() .
  @endreturndesc
@@*/
size_t Util_Strlcat(char* dst, const char* src, size_t dst_size)
{
  const size_t src_len = strlen(src);
  const size_t dst_len = strlen(dst);
  const int dst_remaining = dst_size - dst_len - 1;
  if (dst_remaining > 0)
  {
    strncat(dst, src, dst_remaining);
  }

  return src_len + dst_len;
}

 /*@@
   @routine    Util_StrCmpi
   @date       Mon Jul  5 01:19:00 1999
   @author     Tom Goodale
   @desc
   Case independent strcmp
   @enddesc
   @calls
   @calledby
   @history
   @hdate Wed Oct 13 15:30:57 1999 @hauthor Tom Goodale
   @hdesc Checks the length of the two string first.
   @endhistory
   @var     string1
   @vdesc   First string in comparison
   @vtype   const char *
   @vio     in
   @vcomment

   @endvar
   @var     string2
   @vdesc   Second string in comparison
   @vtype   const char *
   @vio     in
   @vcomment

   @endvar

   @returntype int
   @returndesc
   +ve - string1 > string2
   0   - string1 = string2
   -ve - string1 < string2
   @endreturndesc
@@*/
int Util_StrCmpi (const char *string1, const char *string2)
{
  int retval;


  do
  {
    retval = tolower (*string1) - tolower (*string2);
  } while (! retval && *string1++ && *string2++);

  return (retval);
}

 /*@@
   @routine    Util_StrMemCmpi
   @date       Tue Apr 06 2004
   @author     Erik Schnetter
   @desc
   Case independent strmemcmp: Compare a string against a memory region,
   i.e. a C string against a Fortran string
   @enddesc
   @calls
   @calledby
   @history
   @endhistory
   @var     string1
   @vdesc   First string in comparison (nul-terminated)
   @vtype   const char *
   @vio     in
   @vcomment

   @endvar
   @var     string2
   @vdesc   Second string in comparison (not nul-terminated)
   @vtype   const char *
   @vio     in
   @vcomment

   @endvar
   @var     len2
   @vdesc   Length of the second string
   @vtype   size_t
   @vio     in
   @vcomment

   @endvar

   @returntype int
   @returndesc
   +ve - string1 > string2
   0   - string1 = string2
   -ve - string1 < string2
   @endreturndesc
@@*/
int Util_StrMemCmpi (const char *string1, const char *string2, size_t length2)
{
  int retval;
  size_t last2;


  /* Fortran speciality: ignore trailing blanks */
  last2 = length2;
  while (last2 > 0 && string2[last2-1] == ' ')
  {
    last2--;
  }

  do
  {
    retval = tolower (*string1) - (last2 ? tolower (*string2) : '\0');
  } while (! retval && *string1++ && (string2++, last2--));

  return (retval);
}

 /*@
   @routine    Util_SplitFilename
   @date       Wed Oct 4 10:14:00 2000
   @author     Gabrielle Allen
   @desc
   Splits a filename into its directory and basic filename parts.
   Assigns memory for the two resulting strings, so this should be freed
   when no longer needed.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory
   @var     dir
   @vdesc   The directory part
   @vtype   char **
   @vio     out
   @vcomment

   @endvar
   @var     file
   @vdesc   The file part
   @vtype   char **
   @vio     out
   @vcomment

   @endvar
   @var     string
   @vdesc   The string to split
   @vtype   const char *
   @vio     out
   @vcomment

   @endvar

   @returntype int
   @returndesc
   0  - success
   -1 - out of memory
   @endreturndesc
@@*/
int Util_SplitFilename (char **dir, char **file, const char *string)
{
  char *position;


  *file = Util_Strdup (string);

  if (*file)
  {
    /* Find location of the seperator */
    position = strrchr (*file, '/');
    if (position)
    {
      *dir = *file;
      *position = 0;
      *file = Util_Strdup (position + 1);
    }
    else
    {
      *dir = NULL;
    }
  }

  return (*file ? 0 : -1);
}

 /*@@
   @routine    Util_asprintf
   @date       Thu May 24 16:55:26 2001
   @author     Tom Goodale
   @desc
   Sprintf with memory allocation.  On input
   the buffer should point to a NULL area of memory.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory
   @var     buffer
   @vdesc   Buffer to which to print the string.
   @vtype   char **
   @vio     out
   @vcomment
   *buffer should be NULL on entry.  The routine
   allocates the memory, so the previous contents of
   the pointer are lost.
   On exit the buffer size will be return-value+1 (i.e
   the length of the string plus the \0 ).
   @endvar
   @var     format
   @vdesc   sprintf format string
   @vtype   const char *
   @vio     in
   @vcomment
   This is a standard sprintf format string.
   @endvar
   @var     ...
   @vdesc   Rest of arguments
   @vtype   varargs
   @vio     in
   @vcomment
   These are the arguments necessary for the format string.
   @endvar

   @returntype int
   @returndesc
   The number of bytes written to the buffer.
   @endreturndesc
@@*/
int Util_asprintf(char **buffer, const char *fmt, ...)
{
  int count;
  va_list args;

  va_start(args,fmt);

  count = vsnprintf(NULL, 0, fmt, args);

  va_end(args);

  *buffer = (char *)malloc(count+1);

  if(*buffer)
  {
    va_start(args,fmt);

    vsnprintf(*buffer,count+1,fmt,args);

    va_end(args);
  }
  else
  {
    count = 0;
  }

  return count;
}

 /*@@
   @routine    Util_asprintf
   @date       Thu May 24 16:55:26 2001
   @author     Tom Goodale
   @desc
   Sprintf with memory allocation if necessary.  On input
   the buffer should point to an area of memory of length 'size' .
   @enddesc
   @calls
   @calledby
   @history

   @endhistory
   @var     buffer
   @vdesc   Buffer to which to print the string.
   @vtype   char **
   @vio     out
   @vcomment
   Buffer to which to print string.  If the buffer is too
   small, the buffer is freed and a new buffer big enough to hold
   the string and its null-termination is created.
   @endvar
   @var     size
   @vdesc   initial size of the buffer
   @vtype   int
   @vio     in
   @vcomment
   This is the initial size of the buffer.
   @endvar
   @var     format
   @vdesc   sprintf format string
   @vtype   const char *
   @vio     in
   @vcomment
   This is a standard sprintf format string.
   @endvar
   @var     ...
   @vdesc   Rest of arguments
   @vtype   varargs
   @vio     in
   @vcomment
   These are the arguments necessary for the format string.
   @endvar

   @returntype int
   @returndesc
   The number of bytes written to the buffer.
   @endreturndesc
@@*/
int Util_asnprintf(char **buffer, size_t size, const char *fmt, ...)
{
  size_t count;
  va_list args;

  va_start(args,fmt);

  count = vsnprintf(NULL, 0, fmt, args);

  va_end(args);

  if(count+1 > size)
  {
    /* Use free followed by malloc as realloc may copy memory
     * we are not interested in.
     */
    free(*buffer);
    *buffer = (char *)malloc(count+1);
  }

  if(*buffer)
  {
    va_start(args,fmt);

    vsnprintf(*buffer,count+1,fmt,args);

    va_end(args);
  }
  else
  {
    count = 0;
  }

  return count;
}

int CCTK_CreateDirectory (int mode, const char *pathname)
{
  int retval;
  const char *path;
  char *current;
  const char *token;
  struct stat statbuf;


  current = (char *) malloc (strlen (pathname) + 1);
  if (current)
  {
    retval = 0;
    current[0] = '\0';

    path = pathname;
    while ((token = Util_StrSep (&path, "/")))
    {
      /* Treat first token carefully. */
      if (*current)
      {
        sprintf (current, "%s/%s", current, token);
      }
      else
      {
        strcpy (current, *token ? token : "/");
      }

      if (stat (current, &statbuf))
      {
        if (MKDIR_WRAPPER (current, mode) == -1)
        {
          retval = errno == EEXIST ? 1 : -2;
        }
      }
      else if (! S_ISDIR (statbuf.st_mode))
      {
        retval = -3;
      }
      else
      {
        retval = 1;
      }

      if (retval < 0)
      {
        break;
      }
    }

    if (retval >= 0)
    {
      /* Deal with last component of path */
      if ((size_t) (path - pathname) < strlen (pathname))
      {
        if (stat (pathname, &statbuf))
        {
          retval = 0;
          if (MKDIR_WRAPPER (pathname, mode) == -1)
          {
            retval = errno == EEXIST ? 1 : -2;
          }
        }
        else if (! S_ISDIR (statbuf.st_mode))
        {
          retval = -3;
        }
        else
        {
          retval = 1;
        }
      }
    }

    free (current);

  }
  else
  {
    retval = -1;
  }

  return (retval);
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/



#ifdef TEST_UTIL_STRSEP

#include <stdio.h>

int main(int argc, char *argv[])
{
  const char *argument;
  char *delim;
  const char *token;

  if(argc < 3)
  {
    printf("Usage: %s <string> <delim>\n", argv[0]);
    exit(1);
  }

  argument = argv[1];
  delim = argv[2];

  while((token = Util_StrSep(&argument, delim)))
  {
    printf("Token is     '%s'\n", token);
  }

  if(argument - argv[1] < strlen(argv[1]))
  {
    printf("Remainder is '%s'\n", argument);
  }

  return 0;

}
#endif /*TEST_UTIL_STRSEP */

/******************************************************************************/

#ifdef TEST_UTIL_STRLCPY
/* test_strlcpy -- test driver for Util_Strlcpy() */

#include <string.h>
#include <stdio.h>

/**************************************/

/* prototypes */
size_t tryit(size_t dst_size, const char* src);
void nprint(int n_print, const char* buf);

/* global data structures */
static char buffer[100];

/**************************************/

/*
 * This program is a test driver for  Util_Strlcpy() .
 */

int main(void)
{
size_t n;

n = tryit(9, "hello");
printf("bufsize=9: result=%d buffer=", (int) n);
nprint(9, buffer);

n = tryit(6, "hello");
printf("bufsize=6: result=%d buffer=", (int) n);
nprint(6, buffer);

n = tryit(5, "hello");
printf("bufsize=5: result=%d buffer=", (int) n);
nprint(5, buffer);

n = tryit(4, "hello");
printf("bufsize=4: result=%d buffer=", (int) n);
nprint(4, buffer);

return 0;
}

/**************************************/

size_t tryit(size_t dst_size, const char* src)
{
strcpy(buffer, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");
return Util_Strlcpy(buffer, src, dst_size);
}

/**************************************/

/* print n_print characters of buf[], with visible indication of '\0' */
void nprint(int n_print, const char* buf)
{
  int i;

  printf("\"");
    for (i = 0 ; i < n_print ; ++i)
    {
      if (buf[i] == '\0')
      {
        printf("\\0");
      }
      else
      {
        printf("%c", buf[i]);
      }
    }
  printf("\"\n");
}
#endif /* TEST_UTIL_STRLCPY */

/******************************************************************************/

#ifdef TEST_UTIL_STRLCAT
/* test_strlcat -- test driver for Util_Strlcpy() */

#include <string.h>
#include <stdio.h>

/**************************************/

/* prototypes */
size_t tryit(size_t dst_size, const char* src);
void nprint(int n_print, const char* buf);

/* global data structures */
static char buffer[100];

/**************************************/

/*
 * This program is a test driver for  Util_Strlcpy() .
 */

int main(void)
{
size_t n;

n = tryit(15, "world");
printf("bufsize=15: result=%d buffer=", (int) n);
nprint(20, buffer);

n = tryit(11, "world");
printf("bufsize=11: result=%d buffer=", (int) n);
nprint(20, buffer);

n = tryit(10, "world");
printf("bufsize=10: result=%d buffer=", (int) n);
nprint(20, buffer);

n = tryit(9, "world");
printf("bufsize=9: result=%d buffer=", (int) n);
nprint(20, buffer);

n = tryit(6, "world");
printf("bufsize=6: result=%d buffer=", (int) n);
nprint(20, buffer);

n = tryit(5, "world");
printf("bufsize=5: result=%d buffer=", (int) n);
nprint(20, buffer);

n = tryit(4, "world");
printf("bufsize=4: result=%d buffer=", (int) n);
nprint(20, buffer);

return 0;
}

/**************************************/

size_t tryit(size_t dst_size, const char* src)
{
const char hello[] = "hello\0xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";
memcpy(buffer, hello, sizeof(hello));
return Util_Strlcat(buffer, src, dst_size);
}

/**************************************/

/* print n_print characters of buf[], with visible indication of '\0' */
void nprint(int n_print, const char* buf)
{
  int i;
  int i_null = -1;

  printf("\"");
    for (i = 0 ; i < n_print ; ++i)
    {
      if (buf[i] == '\0')
      {
        if (i_null == -1)
        {
          i_null = i;
        }
        printf("\\0");
      }
      else
      {
        printf("%c", buf[i]);
      }
    }
  if (i_null >= 0)
  {
    printf(" [null at i=%d]", i_null);
  }
  printf("\"\n");
}


#endif /* TEST_UTIL_STRLCAT */

