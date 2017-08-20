/* Table.c -- implementation for key-value tables */
/* $Header$ */

/*@@
 @file          Table.c
 @seeheader     util_Table.h
 @date          Wed Oct 31 16:17:45 MET 2001
 @author        Jonathan Thornburg <jthorn@aei.mpg.de>
 @desc
                This program implements the key-value table API defined
                in util_Table.h, the Cactus Reference Manual, and in
                chapter C of the the Cactus User's Guide.  A slightly
                earlier version of this API is documented in
                  http://www.cactuscode.org/Development/Specs/KeyValueLookup.txt
 @enddesc
 @version       $Id$
 @@*/

/*
 * ***** table of contents for this file *****
 *
 * Growable Array Data Structures
 * Table Data Structures
 * Iterator Data Structures
 * Misc Macros for This File
 * Prototypes for Functions Private to This File
 * Main Table API
 *   Util_TableCreate
 *   Util_TableClone
 *   Util_TableDestroy
 *   Util_TableQueryFlags
 *   Util_TableQueryNKeys
 *   Util_TableQueryMaxKeyLength
 *   Util_TableQueryValueInfo
 *   Util_TableDeleteKey
 *   Util_TableCreateFromString
 *   Util_TableSetFromString
 *   Util_TableSetString
 *   Util_TableGetString
 *   Util_TableSetGeneric
 *   Util_TableSetGenericArray
 *   Util_TableGetGeneric
 *   Util_TableGetGenericArray
 *   Util_TableSet*
 *   Util_TableSet*Array
 *   Util_TableGet*
 *   Util_TableGet*Array
 * Table Iterator API
 *   Util_TableItCreate
 *   Util_TableItClone
 *   Util_TableItDestroy
 *   Util_TableItQueryIsNull
 *   Util_TableItQueryIsNonNull
 *   Util_TableItQueryTableHandle
 *   Util_TableItQueryKeyValueInfo
 *   Util_TableItAdvance
 *   Util_TableItResetToStart
 *   Util_TableItSetToNull
 *   Util_TableItSetToKey
 * Table and Iterator Dump Routines
 *   Util_TablePrintAll
 *   Util_TablePrint
 *   Util_TablePrintPretty
 *   Util_TablePrintAllIterators
 * Internal Support Functions
 *   internal_set
 *   internal_get
 *   get_table_header_ptr
 *   is_bad_key
 *   find_table_entry
 *   insert_table_entry
 *   delete_table_entry_by_key
 *   delete_table_entry_by_ptr
 *   get_iterator_ptr
 *   grow_pointer_array
 *   convert_string_to_number
 */

#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

/* FIXME: C99 defines <stdbool.h>, we should include that or a fake version */
typedef int bool;
#define true    1
#define false   0

#ifndef CCODE
  #define CCODE       /* signal Cactus header files that we're C, not Fortran */
#endif

#include "../AHFD_macros.h"
#include "util_String.h"
#include "util_Table.h"

/*
 * define the following symbol to turn on testing support; this makes
 * the table routines run slower, but may help in catching certain bugs
 */
#undef  UTIL_TABLE_TEST

/* define the following symbol to define Fortran wrappers */
#undef UTIL_TABLE_FORTRAN_WRAPPERS

/* define the following symbol to print various debugging information */
#undef UTIL_TABLE_DEBUG

/* define the following symbol (in addition to UTIL_TABLE_DEBUG) */
/* to print very verbose debugging information */
#undef UTIL_TABLE_DEBUG2

/******************************************************************************/
/***** Growable Array Data Structures *****************************************/
/******************************************************************************/

/*
 * We use "growable arrays" to keep track of all tables and all table
 * iterators.  In both cases we use the same data structure:
 *
 *      int N_objects;          // actual number of tables/iterators
 *      int N_elements;         // actual size of growable array
 *      void *array;            // pointer to malloc-allocated growable array
 *                              // indexed by handle/ihandle
 *
 * Note that the pointer must be  void *  so we can use the  grow_array()
 * function; this pointer should be cast into an actual usable type for
 * normal uses.  Null pointers in the array mark unused array elements.
 */

/*
 * growth policy for growable arrays
 * sequence is
#ifdef UTIL_TABLE_TEST
 *      0, 1, 3, 7, 15, ... entries     (very slow growth
 *                                       ==> better exercise growing code)
#else
 *      0, 10, 30, 70, 150, ... entries
#endif
 * n.b. this grows >= a geometric series
 *      ==> total time in realloc is linear in max array size
 *      (if we just grew in an arithmetic progression then the total
 *       time in realloc() would be quadratic in the max array size)
 */
#ifdef UTIL_TABLE_TEST
  #define GROW(old_n)   (2*(old_n) + 1)
#else
  #define GROW(old_n)   (2*(old_n) + 10)
#endif

/******************************************************************************/
/***** Table Data Structures **************************************************/
/******************************************************************************/

/*
 * The present implementation represents a table as a singly-linked
 * list of table entries.  The code is generally programmed for simplicity,
 * not for maximum performance: linear searches are used everywhere.
 * In practice, we don't expect tables to have very many entries, so
 * this shouldn't be a problem.
 */

struct table_entry
{
  struct table_entry *next;
  char *key;
  int type_code;
  int N_elements;
  void *value;
};

struct table_header
{
  struct table_entry *head;
  int flags;
  int handle;
};

struct scalar_value
{
  int datatype;
  union
  {
    CCTK_INT  int_scalar;
    CCTK_REAL real_scalar;
  } value;
};

/*
 * We keep track of all tables with the following variables
 * (all are static ==> private to this file)
 */

/* number of tables */
static int N_tables = 0;

/* number of elements in the following array */
static int N_thp_array = 0;

/*
 * pointer to growable array of pointers to table headers,
 *            indexed by table handle,
 * with unused array elements set to NULL pointers
 * ... name abbreviates "table-header-pointer array"
 */
void **thp_array = NULL;

/******************************************************************************/
/***** Iterator Data Structures ***********************************************/
/******************************************************************************/

/*
 * This structure represents a table interator.
 *
 * Note that we never modify the table through an iterator,
 * so all the pointers here are to const objects
 */
struct iterator
{
  const struct table_header *thp;   /* must always be non-NULL */
  const struct table_entry *tep;    /* NULL for iterator in */
                                    /* "null-pointer" state */
};

/*
 * We keep track of all iterators with the following variables
 * (all are static ==> private to this file)
 */

/* number of iterators */
static int N_iterators = 0;

/* number of elements in the following array */
static int N_ip_array = 0;

/*
 * pointer to growable array of pointers to iterators,
 *            indexed by iterator handle,
 * with unused array elements set to NULL pointers
 * ... name abbreviates "iterator-pointer array"
 */
void **ip_array = NULL;

/******************************************************************************/
/***** Misc Macros for This File **********************************************/
/******************************************************************************/

#define MIN(x,y)        ((x < y) ? (x) : (y))


/******************************************************************************/
/***** Prototypes for Functions Private to This File **************************/
/******************************************************************************/

/*
 * This is the internal function implementing all the
 *      Util_TableSet*()
 *      Util_TableSet*Array()
 * functions.  It returns their desired return value, i.e.
 *      1 for key was already in table before this call
 *        (old value was replaced)
 *        (it doesn't matter what the old value's type_code and
 *         N_elements were, i.e. these do *not* have to match the
 *         new value),
 *      0 for key was not in table before this call,
 *      UTIL_ERROR_BAD_HANDLE           handle is invalid
 *      UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character
 *      UTIL_ERROR_BAD_INPUT            N_elements < 0
 *      UTIL_ERROR_NO_MEMORY            unable to allocate memory
 */
static
  int internal_set(int handle,
                   int type_code, int N_elements, const void *value,
                   const char *key);

/*
 * This is the internal function implementing all the
 *      Util_TableGet*()
 *      Util_TableGet*Array()
 * functions.  It returns their desired return value, i.e.
 *      number of values stored in  array[]  if ok,
 *      -ve for error, including
 *      UTIL_ERROR_BAD_HANDLE           handle is invalid
 *      UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character
 *      UTIL_ERROR_BAD_INPUT            array != NULL and N_elements < 0
 *      UTIL_ERROR_TABLE_NO_SUCH_KEY    no such key in table
 *      UTIL_ERROR_TABLE_WRONG_DATA_TYPE value has wrong data type
 * If any of the error conditions is returned, the value buffer is unchanged.
 */
static
  int internal_get(int handle,
                   int type_code, int N_value_buffer, void *value_buffer,
                   const char *key);

/* check table handle for validity, return pointer to table header */
static
  struct table_header *get_table_header_ptr(int handle);

/*
 * check if key is syntactically "bad" (eg contains '/' character)
 * returns true for bad key, false for ok
 */
static
  bool is_bad_key(const char *key);

/*
 * delete an entry (specified by its key) from a table
 * return same as Util_TableDeleteKey(), i.e.
 *      0 for ok (key existed before this call, and has now been deleted)
 *      -ve for error, including
 *      UTIL_ERROR_TABLE_NO_SUCH_KEY    no such key in table
 */
static
  int delete_table_entry_by_key(struct table_header *thp, const char *key);

/*
 * delete an entry from a table,
 * specifying the entry by a pointer to the entry *before* it,
 * or NULL to delete the starting entry in the list
 */
static
  void delete_table_entry_by_ptr(struct table_header *thp,
                                 struct table_entry *prev_tep);

/*
 * find table entry for a given key
 * return pointer to it, or NULL if no such key is present in table
 * if  prev_tep_ptr != NULL,
 *         also set *prev_tep_ptr to point to table entry one *before*
 *         the one with the given key, or to NULL if the given key is
 *         the starting entry in the table
 */
static
  struct table_entry *find_table_entry
          (const struct table_header *thp, const char *key,
           struct table_entry **prev_tep_ptr);

/* allocate a new table entry, set its fields to copies of arguments */
static
  int insert_table_entry(struct table_header *thp,
                           const char *key,
                           int type_code, int N_elements, const void *value);

/* check iterator handle for validity, return pointer to iterator */
static
  struct iterator *get_iterator_ptr(int ihandle);

/*
 * This function grows an malloc-allocated array of  void *  pointers
 * via realloc(), initializing the new space to NULL pointers.
 *
 * Arguments:
 * *pN = (in out) array size
 * *pvp_array = (in out) Pointer to growable array of  void *  pointers.
 *
 * Results:
 * This function returns
 *      0 for ok,
 *      -ve for error, including
 *      UTIL_ERROR_NO_MEMORY            can't allocate memory to grow table
 */
static
  int grow_pointer_array(int *pN, void ***pvp_array);

/*
 * This function converts the given string into a scalar value of type
 * CCTK_INT or CCTK_REAL.
 *
 * Arguments:
 * string = (in) null-terminated string to be converted into a number
 * scalar = (out) structure defining the type and value of the number
 */
static
  void convert_string_to_number(const char *string, struct scalar_value *scalar);

/******************************************************************************/
/***** Main Table API *********************************************************/
/******************************************************************************/

/*@@
  @routine      Util_TableCreate
  @desc
                This function creates a new (empty) table.
  @enddesc

  @var          flags
  @vtype        int
  @vdesc        inclusive-or of UTIL_TABLE_FLAGS_* bit flags, must be >= 0
                (n.b. for Fortran users: inclusive-or is the same as sum here,
                since the bit masks are all disjoint)
  @endvar

  @comment
                We require flags >= 0 so other functions can distinguish
                flags from (negative) error codes
  @endcomment

  @returntype   int
  @returndesc
                a handle to the newly-created table,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_NO_MEMORY            unable to allocate memory<BR>
                UTIL_ERROR_TABLE_BAD_FLAGS      flags < 0
  @endreturndesc
  @@*/
int Util_TableCreate(int flags)
{
  #ifdef UTIL_TABLE_DEBUG
  printf("Util_TableCreate()\n");
  #endif

  if (flags < 0)
    return UTIL_ERROR_TABLE_BAD_FLAGS;

  if (N_tables == N_thp_array)
  {
    /* grow  thp_array  to get some room to create the new table */
    #ifdef UTIL_TABLE_DEBUG
    printf("   growing thp_array[] from old size %d\n", N_thp_array);
    #endif
    if (grow_pointer_array(&N_thp_array, &thp_array) < 0)
    {
      return UTIL_ERROR_NO_MEMORY;
    }
    #ifdef UTIL_TABLE_DEBUG
    printf("                         to new size %d\n", N_thp_array);
    #endif
  }

  /* we should now have space to create the new table */
  assert(N_tables < N_thp_array);

  /* find an unused handle */
  #ifdef UTIL_TABLE_DEBUG
  printf("   searching for an unused handle (N_tables=%d N_thp_array=%d)\n",
         N_tables, N_thp_array);
  #endif
    {
  int handle;
  for (handle = 0 ; handle < N_thp_array ; ++handle)
  {
    #ifdef UTIL_TABLE_DEBUG2
    printf("      checking handle=%d\n", handle);
    #endif
    if (thp_array[handle] == NULL)
    {
      /* we've found an unused handle ==> create the table */
      struct table_header *const thp = (struct table_header *)
                                       malloc(sizeof(struct table_header));
      if (thp == NULL)
      {
        return UTIL_ERROR_NO_MEMORY;
      }

      #ifdef UTIL_TABLE_DEBUG
      printf("   using handle=%d\n", handle);
      #endif

      thp->head = NULL;
      thp->flags = flags;
      thp->handle = handle;

      ++N_tables;
      thp_array[handle] = (void *) thp;

      return handle;
    }
  }

  /* we should never get to here! */
  assert(false);
  abort();                                /* internal error (core dump) */
  /* prevent compiler warning 'function should return a value' */
  return(0);
    }
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableCreate)
                          (int *retval, const int *flags);
void CCTK_FCALL CCTK_FNAME(Util_TableCreate)
                          (int *retval, const int *flags)
{
  *retval = Util_TableCreate(*flags);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableClone
  @desc
                This function clones (makes an exact copy of) a table.
                (N.b. the order in which an interator sequences through
                a table may differ in the clone.)
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table to be cloned
  @endvar

  @returntype   int
  @returndesc
                a handle to the clone table, or<BR>
                UTIL_ERROR_NO_MEMORY            unable to allocate memory<BR>
                UTIL_ERROR_TABLE_BAD_FLAGS      flags < 0 in the to-be-cloned
                                                table (this should never happen)
  @endreturndesc
  @@*/
int Util_TableClone(int handle)
{
  const struct table_header *const thp = get_table_header_ptr(handle);
  if (thp == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

    {
  const int clone_handle = Util_TableCreate(thp->flags);
  if (clone_handle < 0)
  {
    return clone_handle;                        /* error creating clone table */
  }

  #ifdef UTIL_TABLE_DEBUG
  printf("Util_TableClone(handle=%d) ==> clone_handle=%d\n",
         handle, clone_handle);
  #endif

  /* copy all the table entries */
    {
  struct table_header *const clone_thp = get_table_header_ptr(clone_handle);
  const struct table_entry *tep;
  for (tep = thp->head ; tep != NULL ; tep = tep->next)
  {
    #ifdef UTIL_TABLE_DEBUG2
    printf("   copying key \"%s\"\n", tep->key);
    #endif
      {
    /* insert_table_entry() does the actual copying */
    int status
        = insert_table_entry(clone_thp,
                             tep->key,
                             tep->type_code, tep->N_elements, tep->value);
    if (status < 0)
    {
      /* this should cleanup as much as we've done so far */
      Util_TableDestroy(clone_handle);
      return status;          /* error inserting table entry into clone table */
    }
      }
  }

  return clone_handle;
    }
    }
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableClone)
                          (int *retval, const int *handle);
void CCTK_FCALL CCTK_FNAME(Util_TableClone)
                          (int *retval, const int *handle)
{
  *retval = Util_TableClone(*handle);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableDestroy
  @desc
                This function destroys a table.
                (Of course, this invalidates any iterators for this table.)
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @returntype   int
  @returndesc
                0 for ok,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid
  @endreturndesc
  @@*/
int Util_TableDestroy(int handle)
{
  struct table_header *const thp = get_table_header_ptr(handle);
  if (thp == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  #ifdef UTIL_TABLE_DEBUG
  printf("Util_TableDestroy(handle=%d)\n", handle);
  #endif

  /* delete all the keys */
  while (thp->head != NULL)
  {
    delete_table_entry_by_ptr(thp, NULL);
  }

  --N_tables;
  thp_array[handle] = NULL;
  free(thp);

  return 0;
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableDestroy)
                          (int *retval, const int *handle);
void CCTK_FCALL CCTK_FNAME(Util_TableDestroy)
                          (int *retval, const int *handle)
{
  *retval = Util_TableDestroy(*handle);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableQueryFlags
  @desc
                This function queries a table's flags word.
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @returntype   int
  @returndesc
                flags if table exists,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid
  @endreturndesc
  @@*/
int Util_TableQueryFlags(int handle)
{
  const struct table_header *const thp = get_table_header_ptr(handle);
  if (thp == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  return thp->flags;
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableQueryFlags)
                          (int *retval, const int *handle);
void CCTK_FCALL CCTK_FNAME(Util_TableQueryFlags)
                          (int *retval, const int *handle)
{
  *retval = Util_TableQueryFlags(*handle);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableQueryNKeys
  @desc
                This function queries the total number of key/value entries
                in a table.
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @returntype   int
  @returndesc
                number of entries (>= 0),<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid
  @endreturndesc
  @@*/
int Util_TableQueryNKeys(int handle)
{
  const struct table_header *const thp = get_table_header_ptr(handle);
  if (thp == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

    {
  int N = 0;
  const struct table_entry *tep;
  for (tep = thp->head ; tep != NULL ; tep = tep->next)
  {
    ++N;
  }

  return N;
    }
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableQueryNKeys)
                          (int *retval, const int *handle);
void CCTK_FCALL CCTK_FNAME(Util_TableQueryNKeys)
                          (int *retval, const int *handle)
{
  *retval = Util_TableQueryNKeys(*handle);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableQueryMaxKeyLength
  @desc
                This function queries the maximum key length in a table.
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @returntype   int
  @returndesc
                maximum key length (>= 0),<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid
  @endreturndesc
  @@*/
int Util_TableQueryMaxKeyLength(int handle)
{
  const struct table_header *const thp = get_table_header_ptr(handle);
  if (thp == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

    {
  int max_length = 0;
  const struct table_entry *tep;
  for (tep = thp->head ; tep != NULL ; tep = tep->next)
  {
    const int length = strlen(tep->key);
    if (length > max_length)
    {
      max_length = length;
    }
  }

  return max_length;
    }
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableQueryMaxKeyLength)
                          (int *retval, const int *handle);
void CCTK_FCALL CCTK_FNAME(Util_TableQueryMaxKeyLength)
                          (int *retval, const int *handle)
{
  *retval = Util_TableQueryMaxKeyLength(*handle);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableQueryValueInfo
  @desc
                This function queries the type and number of elements
                of the value corresponding to a specified key in a table.
                It can also be used to "just" determine whether or not
                a specified key is present in a table.
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          type_code
  @vtype        int *
  @vdesc        pointer to where this function should store
                the value's type code
                (one of the CCTK_VARIABLE_* constants from "cctk_Types.h"),
                or NULL pointer to skip storing this
  @endvar

  @var          N_elements
  @vtype        int *
  @vdesc        pointer to where this function should store
                the number of array elements in the value,
                or NULL pointer to skip storing this
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @returntype   int
  @returndesc
                1 for key is in table,<BR>
                0 for no such key in table
                  (in this case nothing is stored in *type and *N_elements)<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character
  @endreturndesc

  @comment
                Unlike all the other query functions, this function
                returns 0 for no such key in table.  The rationale
                for this design is that by passing NULL pointers for
                type_code and N_elements, this function is then a
                Boolean "is key in table?" predicate.

                If any error code is returned, the user's buffers
                pointed to by type_code and N_elements (if these pointers
                are non-NULL) are unchanged.
  @endcomment
  @@*/
int Util_TableQueryValueInfo(int handle,
                             CCTK_INT *type_code, CCTK_INT *N_elements,
                             const char *key)
{
  const struct table_header *const thp = get_table_header_ptr(handle);
  if (thp == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  if (is_bad_key(key))
  {
    return UTIL_ERROR_TABLE_BAD_KEY;
  }

    {
  const struct table_entry *const tep = find_table_entry(thp, key, NULL);
  if (tep == NULL)
  {
    return 0;                             /* no such key in table */
  }

  if (type_code != NULL)
  {
    *type_code = tep->type_code;
  }
  if (N_elements != NULL)
  {
    *N_elements = tep->N_elements;
  }
  return 1;                               /* key is in table */
    }
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableQueryValueInfo)
                          (int *retval, const int *handle,
                           CCTK_INT *type_code, CCTK_INT *N_elements,
                           ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableQueryValueInfo)
                          (int *retval, const int *handle,
                           CCTK_INT *type_code, CCTK_INT *N_elements,
                           ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableQueryValueInfo(*handle, type_code, N_elements, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableDeleteKey
  @desc
                This function deletes a key (and the corresponding value)
                from a table.

                Note that this invalidates any iterators for this table.
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @returntype   int
  @returndesc
                0 for ok (key existed before this call,
                          and has now been deleted)<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character<BR>
                UTIL_ERROR_TABLE_NO_SUCH_KEY    no such key in table
  @endreturndesc
  @@*/
int Util_TableDeleteKey(int handle, const char *key)
{
  struct table_header *const thp = get_table_header_ptr(handle);
  if (thp == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  if (is_bad_key(key))
  {
    return UTIL_ERROR_TABLE_BAD_KEY;
  }

  return delete_table_entry_by_key(thp, key);
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableDeleteKey)
                          (int *retval, const int *handle, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableDeleteKey)
                          (int *retval, const int *handle, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableDeleteKey(*handle, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableCreateFromString
  @desc
                This function creates a new table (with the case-insensitive
                flag set), and sets values in it based on a string argument.
                The string is interpreted with "parameter-file" semantics.
  @enddesc

  @comment
                The "Implementation Restriction" of Util_TableSetFromString()
                applies here as well.
  @endcomment

  @var          string
  @vtype        const char *
  @vdesc        C-style null-terminated string specifying table contents;
                string has parameter-file semantics
  @endvar

  @returntype   int
  @returndesc
                a handle to the newly-created table,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_NO_MEMORY    unable to allocate memory<BR>
                UTIL_ERROR_BAD_KEY      invalid input: key contains
                                        invalid character<BR>
                UTIL_ERROR_BAD_INPUT    invalid input: can't parse input
                                        string<BR>
                and any error codes returned by
                Util_TableCreate() or Util_TableSetFromString()
  @endreturndesc
  @@*/
int Util_TableCreateFromString(const char string[])
{
  const int handle = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
  if (handle < 0)
  {
    return handle;                                    /* error creating table */
  }

    {
  const int status = Util_TableSetFromString(handle, string);
  if (status < 0)
  {
    Util_TableDestroy(handle);
    return status;                           /* error setting values in table */
  }

  return handle;
    }
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableCreateFromString)
                          (int *retval, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableCreateFromString)
                          (int *retval, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(string)
  *retval = Util_TableCreateFromString(string);
  free(string);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableSetFromString
  @desc
                This function does a sequence of Util_TableSet*() calls
                to set table entries based on a parameter-file--like
                string argument.  For example,
                   Util_TableSetFromString(handle, "order=3 dx=0.1")
                is equivalent to
                   Util_TableSetInt(handle, 3, "order");
                   Util_TableSetReal(handle, 0.1, "dx");
  @enddesc
  @history
  @hdate        Thu 23 May 2002
  @hauthor      Thomas Radke
  @hdesc        Completed for setting string and array values
  @endhistory

  @comment
                Implementation Restriction:<BR>
                The present implementation only recognises integer, real,
                and character-string values (not complex), and integer
                and real arrays.<P>
                In more detail, the strings recognized are defined by the
                following BNF:<BR>
                <BLOCKQUOTE>
                   string -> assign*<BR>
                   assign -> whitespace*<BR>
                   assign -> whitespace* key whitespace* =
                                             whitespace* value delimiter<BR>
                   key    -> any string not containing '/' or '=' or
                             whitespace<BR>
                   value  -> array | int_value | real_value | string_value<BR>
                   array  -> { int_value* } | { real_value* }<BR>
                   int_value    -> anything recognized as a valid integer
                                   by strdol(3) in base 10<BR>
                   real_value   -> anything not recognized as a valid integer
                                   by strtol(3) but recognized as valid by
                                   strdod(3)<BR>
                   string_value -> a C-style string enclosed in "double quotes"
                                   (C-style character escape codes are allowed
                                   ie. '\a', '\b', '\f', '\n', '\r', '\t',
                                       '\v', '\\', '\'', '\"', '\?')<BR>
                   string_value -> A string enclosed in 'single quotes'
                                   (C-style character escape codes are *not*
                                    allowed, ie. every character within the
                                    string is interpreted literally)<BR>
                   delimiter -> end-of-string | whitespace<BR>
                   whitespace --> ' ' | '\t' | '\n' | '\r' | '\f' | '\v'<BR>
                </BLOCKQUOTE>
                where * denotes 0 or more repetitions and | denotes logical or.
                <P>
                Notice also that the keys allowed by this function are
                somewhat more restricted than those allowed by the other
                Util_TableSet*() functions, in that this function disallows
                keys containing '=' and/or whitespace.
  @endcomment

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          string
  @vtype        const char *
  @vdesc        C-style null-terminated string which is parsed as
                described above to determine the keys and values to be
                set in the table.
  @endvar

  @returntype   int
  @returndesc
                the number of successful Util_TableSet*() calls made, or<BR>
                -ve for error, including<BR>
                UTIL_ERROR_NO_MEMORY    unable to allocate memory<BR>
                UTIL_ERROR_BAD_KEY      invalid input: key contains
                                        invalid character<BR>
                UTIL_ERROR_BAD_INPUT    invalid input: can't parse input
                                        string<BR>
                UTIL_ERROR_TABLE_NO_MIXED_TYPE_ARRAY
                                        invalid input: different array elements
                                        differ in their datatypes<BR>
                and any error codes returned by the Util_TableSet*() functions
                Note that in the event of an error return, assignments
                lexicographically earlier in the input string than where
                the error was detected will already have been made in the
                table.  Unfortunately, there is no easy way to find out
                where the error was detected. :(
  @endreturndesc
  @@*/
int CCTK_VarTypeSize (int vtype)
{
  int var_size;


  switch (vtype)
  {
    case CCTK_VARIABLE_BYTE:
      var_size = sizeof (CCTK_BYTE);
      break;

    case CCTK_VARIABLE_INT:
      var_size = sizeof (CCTK_INT);
      break;

    case CCTK_VARIABLE_REAL:
      var_size = sizeof (CCTK_REAL);
      break;


#ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      var_size = sizeof (CCTK_INT1);
      break;
#endif

#ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      var_size = sizeof (CCTK_INT2);
      break;
#endif

#ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      var_size = sizeof (CCTK_INT4);
      break;
#endif

#ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      var_size = sizeof (CCTK_INT8);
      break;
#endif

#ifdef HAVE_CCTK_INT16
    case CCTK_VARIABLE_INT16:
      var_size = sizeof (CCTK_INT16);
      break;
#endif

#ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      var_size = sizeof (CCTK_REAL4);
      break;

    case CCTK_VARIABLE_COMPLEX8:
      var_size = sizeof (CCTK_COMPLEX8);
      break;
#endif

#ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      var_size = sizeof (CCTK_REAL8);
      break;

    case CCTK_VARIABLE_COMPLEX16:
      var_size = sizeof (CCTK_COMPLEX16);
      break;
#endif

#ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      var_size = sizeof (CCTK_REAL16);
      break;

    case CCTK_VARIABLE_COMPLEX32:
      var_size = sizeof (CCTK_COMPLEX32);
      break;
#endif

    case CCTK_VARIABLE_CHAR:
      var_size = sizeof (CCTK_CHAR);
      break;

    case CCTK_VARIABLE_POINTER:
      var_size = sizeof (CCTK_POINTER);
      break;

    case CCTK_VARIABLE_POINTER_TO_CONST:
      var_size = sizeof (CCTK_POINTER_TO_CONST);
      break;

    case CCTK_VARIABLE_FPOINTER:
      var_size = sizeof (CCTK_FPOINTER);
      break;

    default:
      CCTK_CVWarn (4, __LINE__, __FILE__, "Cactus",
                  "CCTK_VarTypeSize: Unknown variable type (%d)", vtype);
      var_size = -1;
  }

  return (var_size);
}


int Util_TableSetFromString(int handle, const char string[])
{
#define WHITESPACE " \t\n\r\f\v"
  struct scalar_value scalar;

  #ifdef UTIL_TABLE_DEBUG
  printf("Util_TableSetFromString(handle=%d, \"%s\")\n", handle, string);
  #endif

  /* make a copy of the string so we can write null characters into it */
  /* to partition it into substrings */
    {
  char *const buffer = Util_Strdup(string);
  if (buffer == NULL)
  {
    return UTIL_ERROR_NO_MEMORY;
  }

    {
  int Set_count = 0, status = 0;
  char *p = buffer;

  while (*p != '\0' && status >= 0)
  {
    /*
     * each pass through this loop processes a single key=value
     * assignment starting at p, creating a table entry for it
     */

    /* skip leading whitespaces */
    p += strspn(p, WHITESPACE);

    #ifdef UTIL_TABLE_DEBUG2
    printf("   skipped over delimiters to p-buffer=%d\n", (int) (p-buffer));
    #endif

    if (*p == '\0')
    {
      break;              /* end of string -- nothing more to do */
    }

      {
    const char *const key = p;                /* key -> "key = value..." */
    char *q = p + strcspn (p, WHITESPACE "=");/* q   -> " = value..." */
    p = q + strspn (q, WHITESPACE);           /* p   -> "= value..." */
    if (*p != '=')
    {
      status = UTIL_ERROR_BAD_INPUT;          /* no '=" in "key=value" string */
      break;
    }

    *q = '\0';                                /* key -> "key" */
    ++p;                                      /* p   -> " value..." */
    p += strspn (p, WHITESPACE);              /* p   -> "value..." */
    if (*p == '\0')
    {
      status = UTIL_ERROR_BAD_INPUT;          /* no value supplied */
      break;
    }

      {
    char *value = p;                          /* value -> "value..." */

    /* split "value..." into "value" and "..." */

    /* check the type of value which is either
     *   - a string enclosed in single or double quotes
     *   - an array of scalars enclosed in round brackets
     *   - a scalar (integer or real)
     */
    if (*value == '\'' || *value == '"' || *value == '{')
    {
      /*
       * this block handles string values and arrays
       */
      q = ((*p == '{') ? "}" : p);            /* q points to delimiter char */

      /* advance to the end of the string or array value */
      do
      {
        /* skip escape character in double-quoted string */
        if (*q == '\"' && *p == '\\' && p[1] != '\0')
        {
          p++;
        }
        p++;
      } while (*p && *p != *q);

      if (*p != *q)
      {
        status = UTIL_ERROR_BAD_INPUT;   /* no closing delimiter found */
        break;                           /* in string or array value */
      }

      /* expand character escape codes in double-quoted string */
      if (*p == '\"')
      {
        while (p > q)
        {
          if (*q == '\\')
          {
            #define CHARACTER_ESCAPE_CODES  "abfnrtv\\\'\"?"
            const char *offset = strchr (CHARACTER_ESCAPE_CODES, q[1]);
            const char character_escape_codes[] =
            {'\a', '\b', '\f', '\n', '\r', '\t', '\v', '\\', '\'', '\"', '\?'};

            if (offset)
            {
              memmove (q, q + 1, p - q);
              q[0] = character_escape_codes[offset - CHARACTER_ESCAPE_CODES];
              q[p - q - 1] = '\0';
            }
            else
            {
              break;                       /* invalid escape code found */
            }
          }
          q++;
        }

        if (p != q)
        {
          status = UTIL_ERROR_BAD_INPUT;   /* invalid escape code found */
          break;
        }
      }
      q = value;            /* q points to the opening delimiter char */
      value++;              /* value skips the delimiter char */
    }
    else
    {
      /*
       * this block handles numbers
       */
      p += strcspn (value, WHITESPACE);
      q = value;
    }

    if (*p != '\0')         /* if we're already at the end of the buffer */
                            /* we don't want to advance further */
    {
      *p++ = '\0';          /* value -> "value", p -> "..." */
    }

    /* set the key/value pair in the table, according to its type */
    if (*q != '{')
    {
      /*
       * this block handles numbers and string values
       */
      if (*q == '\'' || *q == '"')
      {
        status = Util_TableSetString(handle, value, key);
      }
      else
      {
        convert_string_to_number (value, &scalar);
        if (scalar.datatype == CCTK_VARIABLE_INT)
        {
          status = Util_TableSetInt(handle, scalar.value.int_scalar, key);
        }
        else if (scalar.datatype == CCTK_VARIABLE_REAL)
        {
          status = Util_TableSetReal(handle, scalar.value.real_scalar, key);
        }
        else
        {
          status = UTIL_ERROR_BAD_INPUT;   /* can't parse scalar value */
        }
      }

      #ifdef UTIL_TABLE_DEBUG2
      if (status >= 0)
      {
        printf("   ==> storing key='%s', value='%s'\n", key, value);
        printf("   after key=value, advanced p to p-buffer=%d\n",
               (int) (p-buffer));
        printf("   ==> '%s'\n", p);
      }
      #endif
    }
    else
    {
      /*
       * this block handles array values
       */
      int nvals = 0;
      int arraysize = 0;
      int datatype = CCTK_VARIABLE_INT;
      int datatypesize = 0;
      char *array = NULL;


      while (*value)
      {
        /*
         * each pass through this loop processes a single key=value
         * assignment in the array string starting at <value>,
         * extracting the scalar value using convert_string_to_number(),
         * and pushing it to the resulting generic array buffer
         */

        /* skip leading whitespaces */
        value += strspn (value, WHITESPACE);

        /* split "value..." into "value" and "..." */
        q = value + strcspn (value, WHITESPACE);
        if (*q != '\0')         /* if we're already at the end of the list */
                                /* we don't want to advance further */
        {
          *q++ = '\0';          /* value -> "value", q -> "..." */
        }

        convert_string_to_number (value, &scalar);
        if (scalar.datatype == -1)
        {
          datatype = scalar.datatype;
          status = UTIL_ERROR_BAD_INPUT;   /* can't parse array value */
          break;
        }
        if (nvals == 0)
        {
          datatype = scalar.datatype;
          datatypesize = CCTK_VarTypeSize (datatype);
        }
        else if (datatype != scalar.datatype)
        {
          /* all array values must have the same datatype */
          datatype = -1;
          status = UTIL_ERROR_TABLE_NO_MIXED_TYPE_ARRAY;
          break;
        }

        if (nvals >= arraysize)
        {
          if (arraysize == 0)
          {
            /* let's start with an array size of 20 elements */
            arraysize = 20;
            array = malloc (arraysize * datatypesize);
          }
          else
          {
            /* double the array size once it's filled up */
            arraysize *= 2;
            array = realloc (array, arraysize * datatypesize);
          }
          if (array == NULL)
          {
            status = UTIL_ERROR_NO_MEMORY;
            break;
          }
        }

        /* push the new scalar into the array buffer */
        memcpy (array + nvals*datatypesize, &scalar.value, datatypesize);

        #ifdef UTIL_TABLE_DEBUG2
        printf("   ==> storing key='%s', array value='%s'\n", key, value);
        #endif

        nvals++;
        value = q;
      }

      if (datatype == CCTK_VARIABLE_INT || datatype == CCTK_VARIABLE_REAL)
      {
        status = Util_TableSetGenericArray(handle, datatype, nvals, array, key);
      }

      if (array)
      {
        free (array);
      }
    }

    ++Set_count;
      }
      }
  }

  #ifdef UTIL_TABLE_DEBUG2
  printf("   returning with code %d\n", status >= 0 ? Set_count : status);
  #endif

  free(buffer);
  return (status >= 0 ? Set_count : status);
    }
    }
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetFromString)
                          (int *retval, const int *handle, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetFromString)
                          (int *retval, const int *handle, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(string)
  *retval = Util_TableSetFromString(*handle, string);
  free(string);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableSetString
  @desc
                This function sets the value associated with a specified
                key to be (a copy of) a specified character string.

                Note that this invalidates any iterators for this table.
  @enddesc

  @comment
                This function stores the value as array of strlen(string)
                CCTK_CHARs; the stored value does *not* include a terminating
                null character.  (This is convenient for Fortran.)

                The implementation assumes (as is presently the case)
                that a string is in fact an array of CCTK_CHAR, i.e.
                that CCTK_CHAR is the same type as (or at least
                compatible with) char.
  @endcomment

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          string
  @vtype        const char *
  @vdesc        pointer to the (C-style null-terminated) string
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @returntype   int
  @returndesc
                Same as all the other  Util_TableSet*  functions, namely<BR>
                1 for key was already in table before this call
                  (old value was replaced)
                  (it doesn't matter what the old value's type_code and
                   N_elements were, i.e. these do *not* have to match the
                   new value),<BR>
                0 for key was not in table before this call,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character<BR>
                UTIL_ERROR_NO_MEMORY            unable to allocate memory
  @endreturndesc
  @@*/
int Util_TableSetString(int handle,
                        const char *string,
                        const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_CHAR, strlen(string), (const void *) string,
                      key);
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetString)
                          (int *retval, const int *handle, TWO_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetString)
                          (int *retval, const int *handle, TWO_FORTSTRING_ARG)
{
  TWO_FORTSTRING_CREATE(string, key)
  *retval = Util_TableSetString(*handle, string, key);
  free(string); free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableGetString
  @desc
                This function gets a copy of the character-string value
                associated with a specified key, and stores it (or at least
                as much of it as will fit) in a specified character string.
  @enddesc

  @comment
                This function assumes that the value stored in the table
                is an array of CCTK_CHARs, which does *not* include a
                terminating null character.

                The implementation assumes (as is presently the case)
                that a string is in fact an array of CCTK_CHAR, i.e.
                that CCTK_CHAR is the same type as (or at least
                compatible with) char.
  @endcomment

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          buffer_length
  @vtype        int (must be >= 1 if buffer != NULL)
  @vdesc        size of  buffer[]
  @endvar

  @var          buffer
  @vtype        char[]
  @vdesc        a buffer into which this function should store
                (at most  buffer_length-1  characters of) the value,
                terminated by a null character as usual for C strings,
                or NULL pointer to skip storing this
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @returntype   int
  @returndesc
                the string length of the value (as per strlen()),<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character<BR>
                UTIL_ERROR_BAD_INPUT            buffer != NULL
                                                and buffer_length <= 0<BR>
                UTIL_ERROR_TABLE_NO_SUCH_KEY    no such key in table<BR>
                UTIL_ERROR_TABLE_WRONG_DATA_TYPE    value has data type
                                                    other than CCTK_CHAR<BR>
                UTIL_ERROR_TABLE_STRING_TRUNCATED   buffer != NULL and
                                                    value was truncated
                                                    to fit in buffer[]
  @endreturndesc

  @comment
                If the error code UTIL_ERROR_TABLE_STRING_TRUNCATED is
                returned, then the first buffer_length-1 characters of
                the string are returned in the user's buffer (assuming
                buffer is non-NULL), followed by a null character to
                properly terminate the string in the buffer.  If any
                other error code is returned, the user's value buffer
                (pointed to by buffer if this is non-NULL) is unchanged.
  @endcomment
  @@*/
int Util_TableGetString(int handle,
                        int buffer_length, char buffer[],
                        const char *key)
{
  /* string_length = actual length of string, not counting terminating '\0' */
  const int string_length = internal_get(handle,
                                         CCTK_VARIABLE_CHAR,
                                         buffer_length-1, (void *) buffer,
                                         key);
  if (string_length < 0)
  {
    return string_length;                 /* error return from internal_get() */
  }

  /* explicitly add the terminating null character */
  if (buffer != NULL)
  {
    assert(buffer_length >= 1);   /* this should never fail: */
                                  /* internal_get() should return an error */
                                  /* if buffer != NULL and buffer_length <= 0 */
      {
    const int null_posn = MIN(string_length, buffer_length-1);
    buffer[null_posn] = '\0';
      }
  }

  return ((buffer != NULL) && (string_length > buffer_length-1))
         ? UTIL_ERROR_TABLE_STRING_TRUNCATED
         : string_length;
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
/*** FIXME: no fortran wrapper yet ***/
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableSetGeneric
  @desc
                This function sets the value associated with a specified
                key to be a specified value (treated as a 1-element array),
                whose datatype is specified by a CCTK_VARIABLE_* type code.

                Note that this invalidates any iterators for this table.
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          type_code
  @vtype        int
  @vdesc        one of the CCTK_VARIABLE_* constants from "cctk_Types.h",
                describing the actual data type of *value_ptr
  @endvar

  @var          value_ptr
  @vtype        const void *
  @vdesc        a pointer to the value to be associated with the specified key
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @returntype   int
  @returndesc
                1 for key was already in table before this call
                  (old value was replaced)
                  (it doesn't matter what the old value's type_code and
                   N_elements were, i.e. these do *not* have to match the
                   new value),<BR>
                0 for key was not in table before this call,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character<BR>
                UTIL_ERROR_NO_MEMORY            unable to allocate memory
                UTIL_ERROR_BAD_INPUT            unknown  type_code
  @endreturndesc
  @@*/
int Util_TableSetGeneric(int handle,
                         int type_code, const void *value_ptr,
                         const char *key)
{
  return Util_TableSetGenericArray(handle, type_code, 1, value_ptr, key);
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetGeneric)
                          (int *retval, const int *handle,
                           const int *type_code, const CCTK_POINTER *value,
                           ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetGeneric)
                          (int *retval, const int *handle,
                           const int *type_code, const CCTK_POINTER *value,
                           ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetGeneric(*handle, *type_code, value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableSetGenericArray
  @desc
                This function sets the value associated with a specified
                key to be (a copy of) a specified array, whose datatype is
                specified by a CCTK_VARIABLE_* type code.

                Note that this invalidates any iterators for this table.
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          type_code
  @vtype        int
  @vdesc        one of the CCTK_VARIABLE_* constants from "cctk_Types.h",
                describing the actual data type of array[]
  @endvar

  @var          N_elements
  @vtype        int (must be >= 0)
  @vdesc        number of elements in  array[]
  @endvar

  @var          array
  @vtype        const void *
  @vdesc        a pointer to the array (a copy of) which
                is to be associated with the specified key
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @returntype   int
  @returndesc
                1 for key was already in table before this call
                  (old value was replaced)
                  (it doesn't matter what the old value's type_code and
                   N_elements were, i.e. these do *not* have to match the
                   new value),<BR>
                0 for key was not in table before this call,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character<BR>
                UTIL_ERROR_BAD_INPUT            N_elements < 0<BR>
                UTIL_ERROR_NO_MEMORY            unable to allocate memory
  @endreturndesc
  @@*/
int Util_TableSetGenericArray(int handle,
                              int type_code, int N_elements, const void *array,
                              const char *key)
{
  return internal_set(handle,
                      type_code, N_elements, array,
                      key);
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetGenericArray)
                          (int *retval, const int *handle,
                           const int *type_code, const int *N_elements,
                           const CCTK_POINTER array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetGenericArray)
                          (int *retval, const int *handle,
                           const int *type_code, const int *N_elements,
                           const CCTK_POINTER array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetGenericArray(*handle, *type_code, *N_elements,
                                      array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableGetGeneric
  @desc
                This function gets the value of the scalar (1-element array)
                value, or more generally the first array element of the value,
                associated with a specified key.  The value may be of any
                supported datatype; the caller specifies the expected type
                by a CCTK_VARIABLE_* type code.
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          type_code
  @vtype        int
  @vdesc        one of the CCTK_VARIABLE_* constants from "cctk_Types.h",
                describing the expected data type of the table entry.
  @endvar

  @var          value_ptr
  @vtype        void *
  @vdesc        pointer to where this function should store
                a copy of the value associated with the specified key,
                or NULL pointer to skip storing this
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @returntype   int
  @returndesc
                the number of elements in the value,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character<BR>
                UTIL_ERROR_TABLE_NO_SUCH_KEY    no such key in table<BR>
                UTIL_ERROR_TABLE_WRONG_DATA_TYPE value has wrong data type<BR>
                UTIL_ERROR_TABLE_VALUE_IS_EMPTY value is an empty
                                                (0-element) array
  @endreturndesc

  @comment
                Note that it is *not* an error for the value to actually
                be an array with > 1 elements elements; in this case only
                the first element is stored.

                The rationale for this design is that the caller may
                know or suspect that the value is a large array, but
                may only want the first array element; in this case
                this design avoids the caller having to allocate a
                large buffer unnecessarily.

                In contrast, it *is* an error for the value to actually
                be an empty (0-length) array, because then there is no
                ``first array element'' to get.
  @endcomment
  @@*/
int Util_TableGetGeneric(int handle,
                         int type_code, void *value_ptr,
                         const char *key)
{
  const int status
        = Util_TableGetGenericArray(handle, type_code, 1, value_ptr, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableGetGeneric)
                          (int *retval, const int *handle,
                           const int *type_code, CCTK_POINTER *value,
                           ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableGetGeneric)
                          (int *retval, const int *handle,
                           const int *type_code, CCTK_POINTER *value,
                           ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableGetGeneric(*handle, *type_code, value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableGetGenericArray
  @desc
                This is a family of functions, one for each Cactus data type,
                to get a copy of the value associated with a specified key
                (or at least as much of the value as will fit into the
                caller's array).
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          type_code
  @vtype        int
  @vdesc        one of the CCTK_VARIABLE_* constants from "cctk_Types.h",
                describing the expected data type of the table entry.
  @endvar

  @var          N_elements
  @vtype        int (must be >= 0)
  @vdesc        number of elements in  array[]
  @endvar

  @var          array
  @vtype        void *
  @vdesc        a pointer to an array into which this function should store
                (at most  N_elements  elements of) a copy of the value
                associated with the specified key,
                or NULL pointer to skip storing this
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @returntype   int
  @returndesc
                the number of elements in the value,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character<BR>
                UTIL_ERROR_BAD_INPUT            array != NULL and
                                                N_elements < 0<BR>
                UTIL_ERROR_TABLE_NO_SUCH_KEY    no such key in table<BR>
                UTIL_ERROR_TABLE_WRONG_DATA_TYPE value has wrong data type
  @endreturndesc

  @comment
                Note that it is *not* an error for the value to have
                > N_elements elements; in this case only N_elements are
                stored.  The caller can detect this by comparing the
                return value with N_elements.

                The rationale for this design is that the caller may
                know or suspect that the value is a large array, but
                may only want the first few array elements; in this
                case this design avoids the caller having to allocate
                a large buffer unnecessarily.

                It is also *not* an error for the value to have < N_elements
                elements; again the caller can detect this by comparing the
                return value with N_elements.

                Note also that if any error code is returned, the
                caller's value buffer (pointed to by  value_buffer)
                is unchanged.
  @endcomment
  @@*/
int Util_TableGetGenericArray(int handle,
                              int type_code, int N_elements, void *array,
                              const char *key)
{
  return internal_get(handle,
                      type_code, N_elements, array,
                      key);
}

/**************************************/

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableGetGenericArray)
                          (int *retval, const int *handle,
                           const int *type_code, const int *N_elements,
                           CCTK_POINTER array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableGetGenericArray)
                          (int *retval, const int *handle,
                           const int *type_code, const int *N_elements,
                           CCTK_POINTER array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableGetGenericArray(*handle, *type_code, *N_elements,
                                      array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/******************************************************************************/

/*@@
  @routine      Util_TableSet*
  @desc
                This is a family of functions, one for each Cactus data type,
                to set the value associated with a specified key to be a
                specified value (treated as a 1-element array).

                Note that this invalidates any iterators for this table.
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          value
  @vtype        one of
                   CCTK_POINTER, CCTK_POINTER_TO_CONST, CCTK_FPOINTER,
                   CCTK_CHAR,
                   CCTK_BYTE,
                   CCTK_INT, CCTK_INT1, CCTK_INT2, CCTK_INT4, CCTK_INT8,
                   CCTK_INT16,
                   CCTK_REAL, CCTK_REAL4, CCTK_REAL8, CCTK_REAL16,
                   CCTK_COMPLEX, CCTK_COMPLEX8, CCTK_COMPLEX16, CCTK_COMPLEX32
                (not all of these may be supported on any given system)
  @vdesc        the value to be associated with the specified key
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @returntype   int
  @returndesc
                1 for key was already in table before this call
                  (old value was replaced)
                  (it doesn't matter what the old value's type_code and
                   N_elements were, i.e. these do *not* have to match the
                   new value),<BR>
                0 for key was not in table before this call,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character<BR>
                UTIL_ERROR_NO_MEMORY            unable to allocate memory
  @endreturndesc
  @@*/

/**********************************************************/

/*
 * pointers
 */

int Util_TableSetPointer(int handle, CCTK_POINTER value, const char *key)
{
  return Util_TableSetPointerArray(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetPointer)
                          (int *retval, const int *handle,
                           const CCTK_POINTER *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetPointer)
                          (int *retval,
                           const int *handle, const CCTK_POINTER *value,
                           ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetPointer(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

int Util_TableSetPointerToConst(int handle,
                                CCTK_POINTER_TO_CONST value,
                                const char *key)
{
  return Util_TableSetPointerToConstArray(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetPointerToConst)
                          (int *retval, const int *handle,
                           const CCTK_POINTER_TO_CONST *value,
                           ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetPointerTOConst)
                          (int *retval,
                           const int *handle,
                           const CCTK_POINTER_TO_CONST *value,
                           ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetPointerToConst(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

int Util_TableSetFPointer(int handle, CCTK_FPOINTER value, const char *key)
{
  return Util_TableSetFPointerArray(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetFPointer)
                          (int *retval, const int *handle,
                           const CCTK_FPOINTER *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetFPointer)
                          (int *retval, const int *handle,
                           const CCTK_FPOINTER *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetFPointer(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

/*
 * ... the following function (an alias for the previous one) is for
 *     backwards compatability only, and is deprecated as of 4.0beta13
 */
int Util_TableSetFnPointer(int handle, CCTK_FPOINTER value, const char *key)
{
  return Util_TableSetFPointerArray(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetFnPointer)
                          (int *retval, const int *handle,
                           const CCTK_FPOINTER *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetFnPointer)
                          (int *retval, const int *handle,
                           const CCTK_FPOINTER *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetFPointer(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**********************************************************/

/*
 * a single character
 */

int Util_TableSetChar(int handle, CCTK_CHAR value, const char *key)
{
  return Util_TableSetCharArray(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
/*** FIXME: no fortran wrapper yet ***/
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**********************************************************/

/*
 * integers
 */

int Util_TableSetByte(int handle, CCTK_BYTE value, const char *key)
{
  return Util_TableSetByteArray(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetByte)
                          (int *retval, const int *handle,
                           const CCTK_BYTE *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetByte)
                          (int *retval, const int *handle,
                           const CCTK_BYTE *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetByte(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

int Util_TableSetInt(int handle, CCTK_INT value, const char *key)
{
  return Util_TableSetIntArray(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt)
                          (int *retval, const int *handle,
                           const CCTK_INT *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt)
                          (int *retval, const int *handle,
                           const CCTK_INT *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetInt(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

#ifdef HAVE_CCTK_INT1
int Util_TableSetInt1(int handle, CCTK_INT1 value, const char *key)
{
  return Util_TableSetInt1Array(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt1)
                          (int *retval, const int *handle,
                           const CCTK_INT1 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt1)
                          (int *retval, const int *handle,
                           const CCTK_INT1 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetInt1(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif

/**************************************/

#ifdef HAVE_CCTK_INT2
int Util_TableSetInt2(int handle, CCTK_INT2 value, const char *key)
{
  return Util_TableSetInt2Array(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt2)
                          (int *retval, const int *handle,
                           const CCTK_INT2 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt2)
                          (int *retval, const int *handle,
                           const CCTK_INT2 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetInt2(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT2 */

/**************************************/

#ifdef HAVE_CCTK_INT4
int Util_TableSetInt4(int handle, CCTK_INT4 value, const char *key)
{
  return Util_TableSetInt4Array(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt4)
                          (int *retval, const int *handle,
                           const CCTK_INT4 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt4)
                          (int *retval, const int *handle,
                           const CCTK_INT4 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetInt4(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT4 */

/**************************************/

#ifdef HAVE_CCTK_INT8
int Util_TableSetInt8(int handle, CCTK_INT8 value, const char *key)
{
  return Util_TableSetInt8Array(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableSetInt8)
                           (int *retval, const int *handle,
                            const CCTK_INT8 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableSetInt8)
                           (int *retval, const int *handle,
                            const CCTK_INT8 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetInt8(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT8 */

/**************************************/

#ifdef HAVE_CCTK_INT16
int Util_TableSetInt16(int handle, CCTK_INT16 value, const char *key)
{
  return Util_TableSetInt16Array(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableSetInt16)
                           (int *retval, const int *handle,
                            const CCTK_INT16 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableSetInt16)
                           (int *retval, const int *handle,
                            const CCTK_INT16 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetInt16(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT16 */

/**********************************************************/

/*
 * real numbers
 */

int Util_TableSetReal(int handle, CCTK_REAL value, const char *key)
{
  return Util_TableSetRealArray(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetReal)
                          (int *retval, const int *handle,
                           const CCTK_REAL *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetReal)
                          (int *retval, const int *handle,
                           const CCTK_REAL *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetReal(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

#ifdef HAVE_CCTK_REAL4
int Util_TableSetReal4(int handle, CCTK_REAL4 value, const char *key)
{
  return Util_TableSetReal4Array(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetReal4)
                          (int *retval, const int *handle,
                           const CCTK_REAL4 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetReal4)
                          (int *retval, const int *handle,
                           const CCTK_REAL4 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetReal4(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL4 */

/**************************************/

#ifdef HAVE_CCTK_REAL8
int Util_TableSetReal8(int handle, CCTK_REAL8 value, const char *key)
{
  return Util_TableSetReal8Array(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetReal8)
                          (int *retval, const int *handle,
                           const CCTK_REAL8 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetReal8)
                          (int *retval, const int *handle,
                           const CCTK_REAL8 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetReal8(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL8 */

/**************************************/

#ifdef HAVE_CCTK_REAL16
int Util_TableSetReal16(int handle, CCTK_REAL16 value, const char *key)
{
  return Util_TableSetReal16Array(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetReal16)
                          (int *retval, const int *handle,
                           const CCTK_REAL16 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetReal16)
                          (int *retval, const int *handle,
                           const CCTK_REAL16 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetReal16(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL16 */

/**********************************************************/

/*
 * complex numbers
 */


#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplex)
                          (int *retval, const int *handle,
                           const CCTK_COMPLEX *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplex)
                          (int *retval, const int *handle,
                           const CCTK_COMPLEX *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetComplex(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

#ifdef HAVE_CCTK_REAL4
int Util_TableSetComplex8(int handle, CCTK_COMPLEX8 value, const char *key)
{
  return Util_TableSetComplex8Array(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplex8)
                          (int *retval, const int *handle,
                           const CCTK_COMPLEX8 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplex8)
                          (int *retval, const int *handle,
                           const CCTK_COMPLEX8 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetComplex8(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL4 */

/**************************************/

#ifdef HAVE_CCTK_REAL8
int Util_TableSetComplex16(int handle, CCTK_COMPLEX16 value, const char *key)
{
  return Util_TableSetComplex16Array(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplex16)
                          (int *retval, const int *handle,
                           const CCTK_COMPLEX16 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplex16)
                          (int *retval, const int *handle,
                           const CCTK_COMPLEX16 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetComplex16(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL8 */

/**************************************/

#ifdef HAVE_CCTK_REAL16
int Util_TableSetComplex32(int handle, CCTK_COMPLEX32 value, const char *key)
{
  return Util_TableSetComplex32Array(handle, 1, &value, key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplex32)
                          (int *retval, const int *handle,
                           const CCTK_COMPLEX32 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplex32)
                          (int *retval, const int *handle,
                           const CCTK_COMPLEX32 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetComplex32(*handle, *value, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL16 */

/******************************************************************************/

/*@@
  @routine      Util_TableSet*Array
  @desc
                This is a family of functions, one for each Cactus data type,
                to set the value associated with the specified key to be
                (a copy of) a specified array.

                Note that this invalidates any iterators for this table.
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          N_elements
  @vtype        int (must be >= 0)
  @vdesc        number of elements in  array[]
  @endvar

  @var          array
  @vtype        const T[], where T is one of
                   CCTK_POINTER, CCTK_POINTER_TO_CONST, CCTK_FPOINTER,
                   CCTK_CHAR,
                   CCTK_BYTE,
                   CCTK_INT, CCTK_INT1, CCTK_INT2, CCTK_INT4, CCTK_INT8,
                   CCTK_INT16,
                   CCTK_REAL, CCTK_REAL4, CCTK_REAL8, CCTK_REAL16,
                   CCTK_COMPLEX, CCTK_COMPLEX8, CCTK_COMPLEX16, CCTK_COMPLEX32
                (not all of these may be supported on any given system)
  @vdesc        a pointer to the array (a copy of) which
                is to be associated with the specified key
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @returntype   int
  @returndesc
                1 for key was already in table before this call
                  (old value was replaced)
                  (it doesn't matter what the old value's type_code and
                   N_elements were, i.e. these do *not* have to match the
                   new value),<BR>
                0 for key was not in table before this call,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character<BR>
                UTIL_ERROR_BAD_INPUT            N_elements < 0<BR>
                UTIL_ERROR_NO_MEMORY            unable to allocate memory
  @endreturndesc
  @@*/

/**********************************************************/

/*
 * arrays of pointers
 */

int Util_TableSetPointerArray(int handle,
                              int N_elements, const CCTK_POINTER array[],
                              const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_POINTER, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetPointerArray)
                          (int *retval, const int *handle,
                           const int *N_elements, const CCTK_POINTER array[],
                           ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetPointerArray)
                          (int *retval, const int *handle,
                           const int *N_elements, const CCTK_POINTER array[],
                           ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetPointerArray(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

int Util_TableSetPointerToConstArray(int handle,
                                     int N_elements,
                                     const CCTK_POINTER_TO_CONST array[],
                                     const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_POINTER_TO_CONST,
                      N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetPointerToConstArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_POINTER_TO_CONST array[],
                           ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetPointerToConstArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_POINTER_TO_CONST array[],
                           ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetPointerToConstArray(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

int Util_TableSetFPointerArray(int handle,
                               int N_elements, const CCTK_FPOINTER array[],
                               const char *key)
{
  return
    internal_set(handle,
                 CCTK_VARIABLE_FPOINTER, N_elements, (const void *) array,
                 key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetFPointerArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_FPOINTER array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetFPointerArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_FPOINTER array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetFPointerArray(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

/*
 * ... the following function (an alias for the previous one) is for
 *     backwards compatability only, and is deprecated as of 4.0beta13
 */
int Util_TableSetFnPointerArray(int handle,
                                int N_elements, const CCTK_FPOINTER array[],
                                const char *key)
{
  return
    internal_set(handle,
                 CCTK_VARIABLE_FPOINTER, N_elements, (const void *) array,
                 key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetFnPointerArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_FPOINTER array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetFnPointerArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_FPOINTER array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetFPointerArray(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**********************************************************/

/*
 * arrays of characters (i.e. character strings)
 */

int Util_TableSetCharArray(int handle,
                           int N_elements, const CCTK_CHAR array[],
                           const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_CHAR, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetCharArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_CHAR array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetCharArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_CHAR array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetCharArray(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**********************************************************/

/*
 * arrays of integers
 */

int Util_TableSetByteArray(int handle,
                           int N_elements, const CCTK_BYTE array[],
                           const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_BYTE, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetByteArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_BYTE array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetByteArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_BYTE array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetByteArray(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

int Util_TableSetIntArray(int handle,
                          int N_elements, const CCTK_INT array[],
                          const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_INT, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetIntArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_INT array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetIntArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_INT array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetIntArray(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

#ifdef HAVE_CCTK_INT1
int Util_TableSetInt1Array(int handle,
                           int N_elements, const CCTK_INT1 array[],
                           const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_INT1, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableSetInt1Array)
                           (int *retval, const int *handle,
                            const int *N_elements,
                            const CCTK_INT1 array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableSetInt1Array)
                           (int *retval, const int *handle,
                            const int *N_elements,
                            const CCTK_INT1 array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableSetInt1Array (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT1 */

/**************************************/

#ifdef HAVE_CCTK_INT2
int Util_TableSetInt2Array(int handle,
                           int N_elements, const CCTK_INT2 array[],
                           const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_INT2, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt2Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_INT2 array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt2Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_INT2 array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetInt2Array(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT2 */

/**************************************/

#ifdef HAVE_CCTK_INT4
int Util_TableSetInt4Array(int handle,
                           int N_elements, const CCTK_INT4 array[],
                           const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_INT4, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt4Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_INT4 array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt4Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_INT4 array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetInt4Array(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT4 */

/**************************************/

#ifdef HAVE_CCTK_INT8
int Util_TableSetInt8Array(int handle,
                           int N_elements, const CCTK_INT8 array[],
                           const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_INT8, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt8Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_INT8 array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt8Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_INT8 array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetInt8Array(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT8 */

/**************************************/

#ifdef HAVE_CCTK_INT16
int Util_TableSetInt16Array(int handle,
                            int N_elements, const CCTK_INT16 array[],
                            const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_INT16, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt16Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_INT16 array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetInt16Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_INT16 array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetInt16Array(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT16 */

/**********************************************************/

/*
 * arrays of real numbers
 */

int Util_TableSetRealArray(int handle,
                           int N_elements, const CCTK_REAL array[],
                           const char *key)
{
  return
    internal_set(handle,
                    CCTK_VARIABLE_REAL, N_elements, (const void *) array,
                    key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetRealArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_REAL array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetRealArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_REAL array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetRealArray(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

#ifdef HAVE_CCTK_REAL4
int Util_TableSetReal4Array(int handle,
                            int N_elements, const CCTK_REAL4 array[],
                            const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_REAL4, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetReal4Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_REAL4 array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetReal4Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_REAL4 array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetReal4Array(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL4 */

/**************************************/

#ifdef HAVE_CCTK_REAL8
int Util_TableSetReal8Array(int handle,
                            int N_elements, const CCTK_REAL8 array[],
                            const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_REAL8, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetReal8Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_REAL8 array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetReal8Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_REAL8 array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetReal8Array(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL8 */

/**************************************/

#ifdef HAVE_CCTK_REAL16
int Util_TableSetReal16Array(int handle,
                             int N_elements, const CCTK_REAL16 array[],
                             const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_REAL16, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetReal16Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                            const CCTK_REAL16 array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetReal16Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_REAL16 array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetReal16Array(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL16 */

/**********************************************************/

/*
 * arrays of complex numbers
 */

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplexArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_COMPLEX array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplexArray)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_COMPLEX array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetComplexArray(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

#ifdef HAVE_CCTK_REAL4
int Util_TableSetComplex8Array(int handle,
                               int N_elements, const CCTK_COMPLEX8 array[],
                               const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_COMPLEX8, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplex8Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_COMPLEX8 array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableSetComplex8Array)
                           (int *retval, const int *handle,
                            const int *N_elements,
                            const CCTK_COMPLEX8 array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetComplex8Array(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL4 */

/**************************************/

#ifdef HAVE_CCTK_REAL8
int Util_TableSetComplex16Array(int handle,
                                int N_elements, const CCTK_COMPLEX16 array[],
                                const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_COMPLEX16, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplex16Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_COMPLEX16 array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplex16Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_COMPLEX16 array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetComplex16Array(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL8 */

/**************************************/

#ifdef HAVE_CCTK_REAL16
int Util_TableSetComplex32Array(int handle,
                                int N_elements, const CCTK_COMPLEX32 array[],
                                const char *key)
{
  return internal_set(handle,
                      CCTK_VARIABLE_COMPLEX32, N_elements, (const void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplex32Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_COMPLEX32 array[], ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Util_TableSetComplex32Array)
                          (int *retval, const int *handle,
                           const int *N_elements,
                           const CCTK_COMPLEX32 array[], ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key)
  *retval = Util_TableSetComplex32Array(*handle, *N_elements, array, key);
  free(key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL16 */

/******************************************************************************/

/*@@
  @routine      Util_TableGet*
  @desc
                This is a family of functions, one for each Cactus data type,
                to get a copy of the scalar (1-element array) value, or more
                generally the first array element of the value, associated
                with a specified key.
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          value
  @vtype        T *, where T is one of
                   CCTK_POINTER, CCTK_POINTER_TO_CONST, CCTK_FPOINTER,
                   CCTK_CHAR,
                   CCTK_BYTE,
                   CCTK_INT, CCTK_INT1, CCTK_INT2, CCTK_INT4, CCTK_INT8,
                   CCTK_INT16,
                   CCTK_REAL, CCTK_REAL4, CCTK_REAL8, CCTK_REAL16,
                   CCTK_COMPLEX, CCTK_COMPLEX8, CCTK_COMPLEX16, CCTK_COMPLEX32
                (not all of these may be supported on any given system)
  @vdesc        pointer to where this function should store
                a copy of the value associated with the specified key,
                or NULL pointer to skip storing this
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @returntype   int
  @returndesc
                the number of elements in the value,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character<BR>
                UTIL_ERROR_TABLE_NO_SUCH_KEY    no such key in table<BR>
                UTIL_ERROR_TABLE_WRONG_DATA_TYPE value has wrong data type<BR>
                UTIL_ERROR_TABLE_VALUE_IS_EMPTY value is an empty
                                                (0-element) array
  @endreturndesc

  @comment
                Note that it is *not* an error for the value to actually
                be an array with > 1 elements elements; in this case only
                the first element is stored.

                The rationale for this design is that the caller may
                know or suspect that the value is a large array, but
                may only want the first array element; in this case
                this design avoids the caller having to allocate a
                large buffer unnecessarily.

                In contrast, it *is* an error for the value to actually
                be an empty (0-length) array, because then there is no
                ``first array element'' to get.
  @endcomment
  @@*/

/**********************************************************/

/*
 * pointers
 */

int Util_TableGetPointer(int handle, CCTK_POINTER *value, const char *key)
{
  const int status = Util_TableGetPointerArray(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetPointer)
                           (int *retval, const int *handle,
                            CCTK_POINTER *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetPointer)
                           (int *retval, const int *handle,
                            CCTK_POINTER *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetPointer (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

int Util_TableGetPointerToConst(int handle,
                                CCTK_POINTER_TO_CONST *value,
                                const char *key)
{
  const int status = Util_TableGetPointerToConstArray(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetPointerToConst)
                           (int *retval, const int *handle,
                            CCTK_POINTER_TO_CONST *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetPointerToConst)
                           (int *retval, const int *handle,
                            CCTK_POINTER_TO_CONST *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetPointerToConst (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

int Util_TableGetFPointer(int handle, CCTK_FPOINTER *value, const char *key)
{
  const int status = Util_TableGetFPointerArray(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetFPointer)
                           (int *retval, const int *handle,
                            CCTK_FPOINTER *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetFPointer)
                           (int *retval, const int *handle,
                            CCTK_FPOINTER *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetFPointer (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

/*
 * ... the following function (an alias for the previous one) is for
 *     backwards compatability only, and is deprecated as of 4.0beta13
 */
int Util_TableGetFnPointer(int handle, CCTK_FPOINTER *value, const char *key)
{
  const int status = Util_TableGetFPointerArray(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetFnPointer)
                           (int *retval, const int *handle,
                            CCTK_FPOINTER *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetFnPointer)
                           (int *retval, const int *handle,
                            CCTK_FPOINTER *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetFPointer (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**********************************************************/

/*
 * a single character
 */

int Util_TableGetChar(int handle, CCTK_CHAR *value, const char *key)
{
  const int status = Util_TableGetCharArray(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
/*** FIXME: no fortran wrapper yet ***/
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**********************************************************/

/*
 * integers
 */

int Util_TableGetByte(int handle, CCTK_BYTE *value, const char *key)
{
  const int status = Util_TableGetByteArray(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetByte)
                           (int *retval, const int *handle,
                            CCTK_BYTE *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetByte)
                           (int *retval, const int *handle,
                            CCTK_BYTE *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetByte (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

int Util_TableGetInt(int handle, CCTK_INT *value, const char *key)
{
  const int status = Util_TableGetIntArray(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt)
                           (int *retval, const int *handle,
                            CCTK_INT *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt)
                           (int *retval, const int *handle,
                            CCTK_INT *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetInt (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

#ifdef HAVE_CCTK_INT1
int Util_TableGetInt1(int handle, CCTK_INT1 *value, const char *key)
{
  const int status = Util_TableGetInt1Array(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt1)
                           (int *retval, const int *handle,
                            CCTK_INT1 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt1)
                           (int *retval, const int *handle,
                            CCTK_INT1 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetInt1 (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT1 */

/**************************************/

#ifdef HAVE_CCTK_INT2
int Util_TableGetInt2(int handle, CCTK_INT2 *value, const char *key)
{
  const int status = Util_TableGetInt2Array(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt2)
                           (int *retval, const int *handle,
                            CCTK_INT2 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt2)
                           (int *retval, const int *handle,
                            CCTK_INT2 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetInt2 (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT2 */

/**************************************/

#ifdef HAVE_CCTK_INT4
int Util_TableGetInt4(int handle, CCTK_INT4 *value, const char *key)
{
  const int status = Util_TableGetInt4Array(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt4)
                           (int *retval, const int *handle,
                            CCTK_INT4 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt4)
                           (int *retval, const int *handle,
                            CCTK_INT4 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetInt4 (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT4 */

/**************************************/

#ifdef HAVE_CCTK_INT8
int Util_TableGetInt8(int handle, CCTK_INT8 *value, const char *key)
{
  const int status = Util_TableGetInt8Array(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt8)
                           (int *retval, const int *handle,
                            CCTK_INT8 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt8)
                           (int *retval, const int *handle,
                            CCTK_INT8 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetInt8 (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT8 */

/**************************************/

#ifdef HAVE_CCTK_INT16
int Util_TableGetInt16(int handle, CCTK_INT16 *value, const char *key)
{
  const int status = Util_TableGetInt16Array(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt16)
                           (int *retval, const int *handle,
                            CCTK_INT16 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt16)
                           (int *retval, const int *handle,
                            CCTK_INT16 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetInt16 (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT16 */

/**********************************************************/

/*
 * real numbers
 */

int Util_TableGetReal(int handle, CCTK_REAL *value, const char *key)
{
  const int status = Util_TableGetRealArray(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetReal)
                           (int *retval, const int *handle,
                            CCTK_REAL *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetReal)
                           (int *retval, const int *handle,
                            CCTK_REAL *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetReal (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

#ifdef HAVE_CCTK_REAL4
int Util_TableGetReal4(int handle, CCTK_REAL4 *value, const char *key)
{
  const int status = Util_TableGetReal4Array(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetReal4)
                           (int *retval, const int *handle,
                            CCTK_REAL4 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetReal4)
                           (int *retval, const int *handle,
                            CCTK_REAL4 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetReal4 (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL4 */

/**************************************/

#ifdef HAVE_CCTK_REAL8
int Util_TableGetReal8(int handle, CCTK_REAL8 *value, const char *key)
{
  const int status = Util_TableGetReal8Array(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetReal8)
                           (int *retval, const int *handle,
                            CCTK_REAL8 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetReal8)
                           (int *retval, const int *handle,
                            CCTK_REAL8 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetReal8 (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL8 */

/**************************************/

#ifdef HAVE_CCTK_REAL16
int Util_TableGetReal16(int handle, CCTK_REAL16 *value, const char *key)
{
  const int status = Util_TableGetReal16Array(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetReal16)
                           (int *retval, const int *handle,
                            CCTK_REAL16 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetReal16)
                           (int *retval, const int *handle,
                            CCTK_REAL16 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetReal16 (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL16 */

/**********************************************************/

/*
 * complex numbers
 */


#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplex)
                           (int *retval, const int *handle,
                            CCTK_COMPLEX *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplex)
                           (int *retval, const int *handle,
                            CCTK_COMPLEX *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetComplex (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

#ifdef HAVE_CCTK_REAL4
int Util_TableGetComplex8(int handle, CCTK_COMPLEX8 *value, const char *key)
{
  const int status = Util_TableGetComplex8Array(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplex8)
                           (int *retval, const int *handle,
                            CCTK_COMPLEX8 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplex8)
                           (int *retval, const int *handle,
                            CCTK_COMPLEX8 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetComplex8 (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL4 */

/**************************************/

#ifdef HAVE_CCTK_REAL8
int Util_TableGetComplex16(int handle, CCTK_COMPLEX16 *value, const char *key)
{
  const int status = Util_TableGetComplex16Array(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplex16)
                           (int *retval, const int *handle,
                            CCTK_COMPLEX16 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplex16)
                           (int *retval, const int *handle,
                            CCTK_COMPLEX16 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetComplex16 (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL8 */

/**************************************/

#ifdef HAVE_CCTK_REAL16
int Util_TableGetComplex32(int handle, CCTK_COMPLEX32 *value, const char *key)
{
  const int status = Util_TableGetComplex32Array(handle, 1, value, key);
  return (status == 0)
         ? UTIL_ERROR_TABLE_VALUE_IS_EMPTY
         : status;
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplex32)
                           (int *retval, const int *handle,
                            CCTK_COMPLEX32 *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplex32)
                           (int *retval, const int *handle,
                            CCTK_COMPLEX32 *value, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetComplex32 (*handle, value, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL16 */

/******************************************************************************/

/*@@
  @routine      Util_TableGet*Array
  @desc
                This is a family of functions, one for each Cactus data type,
                to get a copy of the value associated with a specified key
                (or at least as much of the value as will fit into the
                caller's array).
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          N_elements
  @vtype        int (must be >= 0)
  @vdesc        number of elements in  array[]
  @endvar

  @var          array
  @vtype        T[], where T is one of
                   CCTK_POINTER, CCTK_POINTER_TO_CONST, CCTK_FPOINTER,
                   CCTK_CHAR,
                   CCTK_BYTE,
                   CCTK_INT, CCTK_INT1, CCTK_INT2, CCTK_INT4, CCTK_INT8,
                   CCTK_INT16,
                   CCTK_REAL, CCTK_REAL4, CCTK_REAL8, CCTK_REAL16,
                   CCTK_COMPLEX, CCTK_COMPLEX8, CCTK_COMPLEX16, CCTK_COMPLEX32
                (not all of these may be supported on any given system)
  @vdesc        an array into which this function should store
                (at most  N_elements  elements of) a copy of the value
                associated with the specified key,
                or NULL pointer to skip storing this
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @comment
                Note that it is *not* an error for the value to have
                > N_elements elements; in this case only N_elements are
                stored.  The caller can detect this by comparing the
                return value with N_elements.

                The rationale for this design is that the caller may
                know or suspect that the value is a large array, but
                may only want the first few array elements; in this
                case this design avoids the caller having to allocate
                a large buffer unnecessarily.

                It is also *not* an error for the value to have < N_elements
                elements; again the caller can detect this by comparing the
                return value with N_elements.

                Note also that if any error code is returned, the
                caller's value buffer (pointed to by  value_buffer)
                is unchanged.
  @endcomment

  @returntype   int
  @returndesc
                the number of elements in the value,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character<BR>
                UTIL_ERROR_BAD_INPUT            array != NULL and N_elements < 0<BR>
                UTIL_ERROR_TABLE_NO_SUCH_KEY    no such key in table<BR>
                UTIL_ERROR_TABLE_WRONG_DATA_TYPE value has wrong data type
  @endreturndesc
  @@*/

/**********************************************************/

/*
 * arrays of pointers
 */

int Util_TableGetPointerArray(int handle,
                              int N_elements, CCTK_POINTER array[],
                              const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_POINTER, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetPointerArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_POINTER array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetPointerArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_POINTER array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetPointerArray (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

int Util_TableGetPointerToConstArray(int handle,
                                     int N_elements,
                                     CCTK_POINTER_TO_CONST array[],
                                     const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_POINTER_TO_CONST,
                      N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetPointerToConstArray)
                           (int *retval, const int *handle,
                            const int *N_elements,
                            CCTK_POINTER_TO_CONST array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetPointerToConstArray)
                           (int *retval, const int *handle,
                            const int *N_elements,
                            CCTK_POINTER_TO_CONST array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetPointerToConstArray(*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

int Util_TableGetFPointerArray(int handle,
                               int N_elements, CCTK_FPOINTER array[],
                               const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_FPOINTER, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetFPointerArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_FPOINTER array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetFPointerArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_FPOINTER array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetFPointerArray (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

/*
 * ... the following function (an alias for the previous one) is for
 *     backwards compatability only, and is deprecated as of 4.0beta13
 */
int Util_TableGetFnPointerArray(int handle,
                                int N_elements, CCTK_FPOINTER array[],
                                const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_FPOINTER, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetFnPointerArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_FPOINTER array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetFnPointerArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_FPOINTER array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetFPointerArray (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**********************************************************/

/*
 * arrays of characters (i.e. character strings)
 */

int Util_TableGetCharArray(int handle,
                           int N_elements, CCTK_CHAR array[],
                           const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_CHAR, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetCharArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_CHAR array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetCharArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_CHAR array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetCharArray (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**********************************************************/

/*
 * arrays of integers
 */

int Util_TableGetByteArray(int handle,
                           int N_elements, CCTK_BYTE array[],
                           const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_BYTE, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetByteArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_BYTE array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetByteArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_BYTE array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetByteArray (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

int Util_TableGetIntArray(int handle,
                          int N_elements, CCTK_INT array[],
                          const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_INT, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetIntArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_INT array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetIntArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_INT array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetIntArray (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

#ifdef HAVE_CCTK_INT1
int Util_TableGetInt1Array(int handle,
                           int N_elements, CCTK_INT1 array[],
                           const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_INT1, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt1Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_INT1 array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt1Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_INT1 array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetInt1Array (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT1 */

/**************************************/

#ifdef HAVE_CCTK_INT2
int Util_TableGetInt2Array(int handle,
                           int N_elements, CCTK_INT2 array[],
                           const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_INT2, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt2Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_INT2 array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt2Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_INT2 array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetInt2Array (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT2 */

/**************************************/

#ifdef HAVE_CCTK_INT4
int Util_TableGetInt4Array(int handle,
                           int N_elements, CCTK_INT4 array[],
                           const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_INT4, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt4Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_INT4 array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt4Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_INT4 array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetInt4Array (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT4 */

/**************************************/

#ifdef HAVE_CCTK_INT8
int Util_TableGetInt8Array(int handle,
                           int N_elements, CCTK_INT8 array[],
                           const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_INT8, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt8Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_INT8 array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt8Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_INT8 array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetInt8Array (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT8 */

/**************************************/

#ifdef HAVE_CCTK_INT16
int Util_TableGetInt16Array(int handle,
                            int N_elements, CCTK_INT16 array[],
                            const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_INT16, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt16Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_INT16 array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetInt16Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_INT16 array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetInt16Array (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_INT16 */

/**********************************************************/

/*
 * arrays of real numbers
 */

int Util_TableGetRealArray(int handle,
                           int N_elements, CCTK_REAL array[],
                           const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_REAL, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetRealArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_REAL array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetRealArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_REAL array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetRealArray (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

#ifdef HAVE_CCTK_REAL4
int Util_TableGetReal4Array(int handle,
                            int N_elements, CCTK_REAL4 array[],
                            const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_REAL4, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetReal4Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_REAL4 array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetReal4Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_REAL4 array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetReal4Array (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL4 */

/**************************************/

#ifdef HAVE_CCTK_REAL8
int Util_TableGetReal8Array(int handle,
                            int N_elements, CCTK_REAL8 array[],
                            const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_REAL8, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetReal8Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_REAL8 array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetReal8Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_REAL8 array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetReal8Array (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL8 */

/**************************************/

#ifdef HAVE_CCTK_REAL16
int Util_TableGetReal16Array(int handle,
                             int N_elements, CCTK_REAL16 array[],
                             const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_REAL16, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetReal16Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_REAL16 array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetReal16Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_REAL16 array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetReal16Array (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL16 */

/**********************************************************/


#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplexArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_COMPLEX array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplexArray)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_COMPLEX array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetComplexArray (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */

/**************************************/

#ifdef HAVE_CCTK_REAL4
int Util_TableGetComplex8Array(int handle,
                               int N_elements, CCTK_COMPLEX8 array[],
                               const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_COMPLEX8, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplex8Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_COMPLEX8 array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplex8Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_COMPLEX8 array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetComplex8Array (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL4 */

/**************************************/

#ifdef HAVE_CCTK_REAL8
int Util_TableGetComplex16Array(int handle,
                                int N_elements, CCTK_COMPLEX16 array[],
                                const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_COMPLEX16, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplex16Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_COMPLEX16 array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplex16Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_COMPLEX16 array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetComplex16Array (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL8 */

/**************************************/

#ifdef HAVE_CCTK_REAL16
int Util_TableGetComplex32Array(int handle,
                                int N_elements, CCTK_COMPLEX32 array[],
                                const char *key)
{
  return internal_get(handle,
                      CCTK_VARIABLE_COMPLEX32, N_elements, (void *) array,
                      key);
}

#ifdef UTIL_TABLE_FORTRAN_WRAPPERS
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplex32Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_COMPLEX32 array[],
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (Util_TableGetComplex32Array)
                           (int *retval, const int *handle,
                            const int *N_elements, CCTK_COMPLEX32 array[],
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (key)
  *retval = Util_TableGetComplex32Array (*handle, *N_elements, array, key);
  free (key);
}
#endif  /* UTIL_TABLE_FORTRAN_WRAPPERS */
#endif  /* CCTK_REAL16 */

/******************************************************************************/
/***** Table Iterator API *****************************************************/
/******************************************************************************/

/*@@
  @routine      Util_TableItCreate
  @desc
                This function creates a new table iterator.  The iterator
                points to the starting entry in the table's traversal order.
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @returntype   int
  @returndesc
                a handle to the newly-created iterator,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           table handle is invalid<BR>
                UTIL_ERROR_NO_MEMORY            unable to allocate memory
  @endreturndesc
  @@*/
int Util_TableItCreate(int handle)
{
  const struct table_header *const thp = get_table_header_ptr(handle);
  if (thp == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  #ifdef UTIL_TABLE_DEBUG
  printf("Util_TableItCreate(handle=%d)\n", handle);
  #endif

  if (N_iterators == N_ip_array)
  {
    /* grow  iterator_array  to get some room to create the new table */
    #ifdef UTIL_TABLE_DEBUG
    printf("   growing ip_array[] from old size %d\n",
           N_ip_array);
    #endif
    if (grow_pointer_array(&N_ip_array, &ip_array) < 0)
    {
      return UTIL_ERROR_NO_MEMORY;                        /* can't grow array */
    }
    #ifdef UTIL_TABLE_DEBUG
    printf("      to new size %d\n",
           N_ip_array);
    #endif
  }

  /* we should now have space to create the new iterator */
  assert(N_iterators < N_ip_array);

  /* find an unused iterator handle */
  #ifdef UTIL_TABLE_DEBUG
  printf("   searching for an unused iterator handle\n");
  printf("   (N_iterators=%d N_ip_array=%d\n", N_iterators, N_ip_array);
  #endif
    {
  int ihandle;
    for (ihandle = 0 ; ihandle < N_ip_array ; ++ihandle)
    {
      #ifdef UTIL_TABLE_DEBUG2
      printf("      checking ihandle=%d\n", ihandle);
      #endif
      if (ip_array[ihandle] == NULL)
      {
        /* we've found an unused ihandle ==> create the iterator */
        struct iterator *const ip = (struct iterator *)
                                    malloc(sizeof(struct iterator));
        if (ip == NULL)
        {
          return UTIL_ERROR_NO_MEMORY;         /* can't allocate new iterator */
        }

        #ifdef UTIL_TABLE_DEBUG2
        printf("   using ihandle=%d\n", ihandle);
        #endif

        ip->thp = thp;
        ip->tep = thp->head;    /* iterator initially -> start of table */

        ++N_iterators;
        ip_array[ihandle] = (void *) ip;

        return ihandle;                                      /* NORMAL RETURN */
      }
    }

  /* we should never get to here! */
  assert(false);
  abort();                                      /* internal error (core dump) */
  /* prevent compiler warning 'function should return a value' */
  return(0);
    }
}

/******************************************************************************/

/*@@
  @routine      Util_TableItClone
  @desc
                This function clones (makes an exact copy of) a table
                iterator.  That is, it creates a new iterator which points
                to the same table entry as an existing iterator.
  @enddesc

  @var          ihandle
  @vtype        int
  @vdesc        handle to the iterator to be cloned
  @endvar

  @returntype   int
  @returndesc
                a handle to the newly-created iterator,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           iterator handle is invalid<BR>
                UTIL_ERROR_NO_MEMORY            unable to allocate memory
  @endreturndesc
  @@*/
int Util_TableItClone(int ihandle)
{
  struct iterator *const ip = get_iterator_ptr(ihandle);
  if (ip == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

    {
  const int clone_ihandle = Util_TableItCreate(ip->thp->handle);
  if (clone_ihandle < 0)
  {
    return clone_ihandle;                 /* error in creating clone iterator */
  }

    {
  struct iterator *const clone_ip = get_iterator_ptr(clone_ihandle);
  clone_ip->tep = ip->tep;
  return clone_ihandle;
    }
    }
}

/******************************************************************************/

/*@@
  @routine      Util_TableItDestroy
  @desc
                This function destroys a table iterator.
  @enddesc

  @var          ihandle
  @vtype        int
  @vdesc        handle to the iterator
  @endvar

  @returntype   int
  @returndesc
                0 for ok,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           iterator handle is invalid
  @endreturndesc
  @@*/
int Util_TableItDestroy(int ihandle)
{
  struct iterator *const ip = get_iterator_ptr(ihandle);
  if (ip == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  #ifdef UTIL_TABLE_DEBUG
  printf("Util_TableItDestroy(ihandle=%d)\n", ihandle);
  #endif

  --N_iterators;
  ip_array[ihandle] = NULL;
  free(ip);

  return 0;                               /* ok */
}

/******************************************************************************/

/*@@
  @routine      Util_TableItQueryIsNull
  @desc
                This function queries whether a table iterator is in the
                "null-pointer" state, i.e. whether it does *not* point
                to some table entry.

                Bad things (garbage results, core dumps) may happen if
                you call this function on a table iterator which has been
                invalidated by a change in the table's contents.
  @enddesc

  @var          ihandle
  @vtype        int
  @vdesc        handle to the iterator
  @endvar

  @returntype   int
  @returndesc
                1 for iterator is in "null-pointer" state,<BR>
                0 for iterator points to some table entry,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           iterator handle is invalid
  @endreturndesc
  @@*/
int Util_TableItQueryIsNull(int ihandle)
{
  const struct iterator *const ip = get_iterator_ptr(ihandle);
  if (ip == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  return (ip->tep == NULL)
         ? 1                              /* iterator in "null-pointer" state */
         : 0;                             /* iterator -> some table entry */
}

/******************************************************************************/

/*@@
  @routine      Util_TableItQueryIsNonNull
  @desc
                This function queries whether a table iterator is *not* in
                the "null-pointer" state, i.e. whether it points to some
                table entry.

                Bad things (garbage results, core dumps) may happen if
                you call this function on an iterator which has been
                invalidated by a change in the table's contents.
  @enddesc

  @var          ihandle
  @vtype        int
  @vdesc        handle to the iterator
  @endvar

  @returntype   int
  @returndesc
                1 for iterator points to some table entry,<BR>
                0 for iterator is in "null-pointer" state,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           iterator handle is invalid
  @endreturndesc
  @@*/
int Util_TableItQueryIsNonNull(int ihandle)
{
  const struct iterator *const ip = get_iterator_ptr(ihandle);
  if (ip == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  return (ip->tep == NULL)
         ? 0                              /* iterator in "null-pointer" state */
         : 1;                             /* iterator -> some table entry */
}

/******************************************************************************/

/*@@
  @routine      Util_TableItQueryTableHandle
  @desc
                This function queries which table a table iterator points
                into.

                Note that this is always well-defined, even if the iterator
                is in the "null-pointer" state, and even if the iterator
                has been invalidated by a change in the table's contents.
  @enddesc

  @var          ihandle
  @vtype        int
  @vdesc        handle to the iterator
  @endvar

  @returntype   int
  @returndesc
                table handle,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           iterator handle is invalid
  @endreturndesc
  @@*/
int Util_TableItQueryTableHandle(int ihandle)
{
  const struct iterator *const ip = get_iterator_ptr(ihandle);
  if (ip == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  return ip->thp->handle;
}

/******************************************************************************/

/*@@
  @routine      Util_TableItQueryKeyValueInfo
  @desc
                This function queries the key and the type and number of
                elements of the value corresponding to that key, of the
                table entry to which an iterator points.  This is in fact
                the main purpose of iterators.

                Bad things (garbage results, core dumps) may happen if
                you call this function on an iterator which has been
                invalidated by a change in the table's contents.
  @enddesc

  @var          ihandle
  @vtype        int
  @vdesc        handle to the iterator
  @endvar

  @var          key_buffer_length,
  @vtype        int (must be >= 1 if key_buffer != NULL)
  @vdesc        length of  key_buffer[]  buffer
  @endvar

  @var          key_buffer,
  @vtype        char []
  @vdesc        a buffer into which this function should store
                (at most  key_buffer_length-1  characters of) the key,
                terminated by a null character as usual for C strings,
                or NULL pointer to skip storing this
  @endvar

  @var          type_code
  @vtype        CCTK_INT *
  @vdesc        pointer to where this function should store
                the value's type code
                (one of the CCTK_VARIABLE_* constants from "cctk_Types.h"),
                or NULL pointer to skip storing this
  @endvar

  @var          N_elements
  @vtype        CCTK_INT *
  @vdesc        pointer to where this function should store
                the number of array elements in the value,
                or NULL pointer to skip storing this
  @endvar

  @returntype   int
  @returndesc
                the string length of the key (as per strlen()),<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           iterator handle is invalid<BR>
                UTIL_ERROR_TABLE_ITERATOR_IS_NULL  iterator is in
                                                   "null-pointer" state<BR>
                UTIL_ERROR_TABLE_STRING_TRUNCATED  key_buffer != NULL and
                                                   key was truncated
                                                   to fit in key_buffer[]
  @endreturndesc

  @comment
                If the error code UTIL_ERROR_TABLE_STRING_TRUNCATED is
                returned, then the first key_buffer_length-1 characters of
                the string are returned in the user's key buffer (assuming
                key_buffer is non-NULL), followed by a null character to
                properly terminate the string in the buffer.  If any
                other error code is returned, the user's key buffer
                (pointed to by key_buffer if this is non-NULL) is unchanged.
  @endcomment
  @@*/
int Util_TableItQueryKeyValueInfo(int ihandle,
                                  int key_buffer_length, char key_buffer[],
                                  CCTK_INT *type_code, CCTK_INT *N_elements)
{
  const struct iterator *const ip = get_iterator_ptr(ihandle);
  if (ip == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

    {
  const struct table_entry *const tep = ip->tep;
  if (tep == NULL)
  {
    return UTIL_ERROR_TABLE_ITERATOR_IS_NULL;
  }

    {
  const int actual_key_length = strlen(tep->key);

  /* store the fixed-length output arguments first, so the caller */
  /* will have them even if we hit an error trying to copy the key */
  if (type_code != NULL)
  {
    *type_code = tep->type_code;
  }
  if (N_elements != NULL)
  {
    *N_elements = tep->N_elements;
  }

  if (key_buffer != NULL)
  {
    const int N_key_copy = MIN(key_buffer_length-1, actual_key_length);
    if (N_key_copy < 0)     /* can only happen if key_buffer_length <= 0 */
    {
      /*
       * We have to bail out now, before trying the memcpy(), because
       * memcpy() takes a size_t (= unsigned) value for its count of how
       * many chars to copy, and converting our -ve N_key_copy to size_t
       * would give a huge +ve count :( :(
       */
      return UTIL_ERROR_TABLE_STRING_TRUNCATED;
    }
    memcpy(key_buffer, tep->key, N_key_copy);
    key_buffer[N_key_copy] = '\0';
    if (N_key_copy < actual_key_length)
    {
      return UTIL_ERROR_TABLE_STRING_TRUNCATED;
    }
  }

  return actual_key_length;         /* ok */
    }
    }
}

/******************************************************************************/

/*@@
  @routine      Util_TableItAdvance
  @desc
                This function advances a table iterator to the next entry
                in the table's traversal order.

                Bad things (garbage results, core dumps) may happen if
                you call this function on an iterator which has been
                invalidated by a change in the table's contents.
  @enddesc

  @var          ihandle
  @vtype        int
  @vdesc        handle to the iterator
  @endvar

  @returntype   int
  @returndesc
                same as that of Util_TableItQueryNonNull(ihandle)
                after advancing the iterator, i.e.<BR>
                1 for ok and iterator now points to some table element,<BR>
                0 for advance-past-last-entry
                  (sets iterator to "null-pointer" state),<BR>
                0 if iterator was already in "null-pointer" state)
                  (in this case this call is a no-op),<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           iterator handle is invalid
  @endreturndesc
  @@*/
int Util_TableItAdvance(int ihandle)
{
  struct iterator *const ip = get_iterator_ptr(ihandle);
  if (ip == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  if (ip->tep == NULL)
  {
    return 0;         /* iterator was already in "null-pointer" state */
  }

  ip->tep = ip->tep->next;

  return (ip->tep == NULL)
         ? 0              /* advance past last entry */
                          /* ==> iterator now in "null-pointer" state */
         : 1;             /* ok */
}

/******************************************************************************/

/*@@
  @routine      Util_TableItResetToStart
  @desc
                This function resets a table iterator to point to the
                starting entry in the table's traversal order.

                Note that it is always ok to call this function, even
                if the iterator has been invalidated by a change in the
                table's contents.
  @enddesc

  @var          ihandle
  @vtype        int
  @vdesc        handle to the iterator
  @endvar

  @returntype   int
  @returndesc
                same as that of Util_TableItQueryNonNull(ihandle)
                after resetting the iterator, i.e.<BR>
                1 for ok and iterator now points to some table element,<BR>
                0 for ok and iterator is now in "null-pointer" state
                  (means table is empty)<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           iterator handle is invalid
  @endreturndesc
  @@*/
int Util_TableItResetToStart(int ihandle)
{
  struct iterator *const ip = get_iterator_ptr(ihandle);
  if (ip == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  ip->tep = ip->thp->head;
  return (ip->tep == NULL)
         ? 0              /* ok, iterator is now in "null-pointer" state */
                          /*     (table must be empty) */
         : 1;             /* ok, iterator points to some table element */
}

/******************************************************************************/

/*@@
  @routine      Util_TableItSetToNull
  @desc
                This function sets a table iterator to the "null-pointer"
                state.

                Note that it is always ok to call this function, even
                if the iterator has been invalidated by a change in the
                table's contents.
  @enddesc

  @var          ihandle
  @vtype        int
  @vdesc        handle to the iterator
  @endvar

  @returntype   int
  @returndesc
                0 for ok,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           iterator handle is invalid
  @endreturndesc
  @@*/
int Util_TableItSetToNull(int ihandle)
{
  struct iterator *const ip = get_iterator_ptr(ihandle);
  if (ip == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  ip->tep = NULL;
  return 0;                               /* ok */
}

/******************************************************************************/

/*@@
  @routine      Util_TableItSetToKey
  @desc
                This function sets a table iterator to point to a
                specified table entry.  It has the same effect as
                Util_TableItResetToStart() followed by repeatedly
                calling Util_TableItAdvance() until the iterator
                points to the desired table entry.

                Note that it is always ok to call this function, even
                if the iterator has been invalidated by a change in the
                table's contents.
  @enddesc

  @var          ihandle
  @vtype        int
  @vdesc        handle to the iterator
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @returntype   int
  @returndesc
                0 for ok,<BR>
                UTIL_ERROR_BAD_HANDLE           iterator handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character<BR>
                UTIL_ERROR_TABLE_NO_SUCH_KEY    no such key in table
  @endreturndesc
  @@*/
int Util_TableItSetToKey(int ihandle, const char *key)
{
  struct iterator *const ip = get_iterator_ptr(ihandle);
  if (ip == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  if (is_bad_key(key))
  {
    return UTIL_ERROR_TABLE_BAD_KEY;
  }

  ip->tep = find_table_entry(ip->thp, key, NULL);
  if (ip->tep == NULL)
  {
    return UTIL_ERROR_TABLE_NO_SUCH_KEY;
  }

  return 0;
}

/******************************************************************************/
/***** Internal Support Functions *********************************************/
/******************************************************************************/

/*@@
  @routine      internal_set
  @desc
                This is the internal function implementing all the
                        Util_TableSet*()
                        Util_TableSet*Array()
                functions except Util_TableSetString().  It sets the
                value associated with a specified key, to be a copy
                of a specified array.

                Note that this invalidates any iterators for this table.
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          N_elements
  @vtype        int (must be >= 0)
  @vdesc        number of elements in  array[]
  @endvar

  @var          array
  @vtype        const T[], where T is one of
                   CCTK_POINTER, CCTK_POINTER_TO_CONST, CCTK_FPOINTER,
                   CCTK_CHAR,
                   CCTK_BYTE,
                   CCTK_INT, CCTK_INT1, CCTK_INT2, CCTK_INT4, CCTK_INT8,
                   CCTK_INT16,
                   CCTK_REAL, CCTK_REAL4, CCTK_REAL8, CCTK_REAL16,
                   CCTK_COMPLEX, CCTK_COMPLEX8, CCTK_COMPLEX16, CCTK_COMPLEX32
                (not all of these may be supported on any given system)
  @vdesc        the array (a copy of) which is to be
                associated with the specified key
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @returntype   int
  @returndesc
                1 for key was already in table before this call
                  (old value was replaced)
                  (it doesn't matter what the old value's type_code and
                   N_elements were, i.e. these do *not* have to match the
                   new value),<BR>
                0 for key was not in table before this call,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character<BR>
                UTIL_ERROR_BAD_INPUT            N_elements < 0<BR>
                                                or unknown  type_code
                UTIL_ERROR_NO_MEMORY            unable to allocate memory
  @endreturndesc
  @@*/
static
  int internal_set(int handle,
                   int type_code, int N_elements, const void *value,
                   const char *key)
{
  #ifdef UTIL_TABLE_DEBUG
  printf("internal_set(handle=%d, type_code=%d, N_elements=%d, key=\"%s\")\n",
         handle, type_code, N_elements, key);
  #endif

    {
  struct table_header *const thp = get_table_header_ptr(handle);
  if (thp == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  if (is_bad_key(key))
  {
    return UTIL_ERROR_TABLE_BAD_KEY;
  }
  if (N_elements < 0)
  {
    return UTIL_ERROR_BAD_INPUT;
  }

  /* if key is already in table, delete it */
  /* ... this is a harmless no-op if it's not already in the table */
    {
  int return_value;
  switch (delete_table_entry_by_key(thp, key))
  {
    case 0:
      return_value = 1;       /* key was already in table before this call */
                              /* (we've just deleted it, and we're about */
                              /*  to set the replacement in the table) */
      break;
    case UTIL_ERROR_TABLE_NO_SUCH_KEY:
      return_value = 0;       /* key was not in table before this call */
      break;
    default:
      /* unexpected return code from  delete_table_entry_by_key() */
      /* (this should never happen!) */
      assert(false);
      abort();                                  /* internal error (core dump) */
  }

    {
  const int status = insert_table_entry(thp,
                                        key,
                                        type_code, N_elements, value);
  if (status < 0)
  {
    return status;                        /* error inserting entry into table */
  }

  return return_value;
    }
    }
    }
}

/******************************************************************************/

/*@@
  @routine      internal_get
  @desc
                This is the internal function implementing all the
                        Util_TableGet*()
                        Util_TableGet*Array()
                functions except for Util_TableGetString().  It copies
                up to N_elements of the value associated with a specified
                key, into a user-supplied buffer.
  @enddesc

  @var          handle
  @vtype        int
  @vdesc        handle to the table
  @endvar

  @var          N_value_buffer
  @vtype        int (must be >= 0)
  @vdesc        number of elements in  array[]
  @endvar

  @var          value_buffer
  @vtype        T[], where T is one of
                   CCTK_POINTER, CCTK_POINTER_TO_CONST, CCTK_FPOINTER,
                   CCTK_CHAR,
                   CCTK_BYTE,
                   CCTK_INT, CCTK_INT1, CCTK_INT2, CCTK_INT4, CCTK_INT8,
                   CCTK_INT16,
                   CCTK_REAL, CCTK_REAL4, CCTK_REAL8, CCTK_REAL16,
                   CCTK_COMPLEX, CCTK_COMPLEX8, CCTK_COMPLEX16, CCTK_COMPLEX32
                (not all of these may be supported on any given system)
  @vdesc        an array into which this function should store
                (at most  N_elements  elements of) a copy of the value
                associated with the specified key,
                or NULL pointer to skip storing this
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to the key (a C-style null-terminated string)
  @endvar

  @returntype   int
  @returndesc
                number of elements in the value,<BR>
                -ve for error, including<BR>
                UTIL_ERROR_BAD_HANDLE           handle is invalid<BR>
                UTIL_ERROR_TABLE_BAD_KEY        key contains '/' character<BR>
                UTIL_ERROR_BAD_INPUT            N_value_buffer < 0<BR>
                UTIL_ERROR_BAD_INPUT            value_buffer != NULL
                                                and N_value_buffer < 0<BR>
                UTIL_ERROR_TABLE_NO_SUCH_KEY    no such key in table<BR>
                UTIL_ERROR_TABLE_WRONG_DATA_TYPE value has wrong data type
  @endreturndesc

  @comment
                Note that it is *not* an error for the value to have
                > N_value_buffer elements; in this case only N_value_buffer
                are stored.  The caller can detect this by comparing the
                return value with N_value_buffer.

                Note also that if any error code is returned, the
                caller's value buffer (pointed to by  value_buffer)
                is unchanged.
  @endcomment
  @@*/
static
  int internal_get(int handle,
                   int type_code, int N_value_buffer, void *value_buffer,
                   const char *key)
{
  #ifdef UTIL_TABLE_DEBUG
  printf(
     "internal_get(handle=%d, type_code=%d, N_value_buffer=%d, key=\"%s\")\n",
         handle, type_code, N_value_buffer, key);
  #endif

    {
  const struct table_header *const thp = get_table_header_ptr(handle);
  if (thp == NULL)
  {
    return UTIL_ERROR_BAD_HANDLE;
  }

  if (is_bad_key(key))
  {
    return UTIL_ERROR_TABLE_BAD_KEY;
  }

    {
  const struct table_entry *const tep = find_table_entry(thp, key, NULL);
  if (tep == NULL)
  {
    return UTIL_ERROR_TABLE_NO_SUCH_KEY;              /* no such key in table */
  }

  if (tep->type_code != type_code)
  {
    return UTIL_ERROR_TABLE_WRONG_DATA_TYPE;     /* value has wrong data type */
  }

  if (value_buffer != NULL)
  {
    if (N_value_buffer < 0)
    {
      return UTIL_ERROR_BAD_INPUT;
    }
      {
    const int N_copy = MIN(N_value_buffer, tep->N_elements);
    const size_t sizeof_N_copy_elements = N_copy * CCTK_VarTypeSize(type_code);
    #ifdef UTIL_TABLE_DEBUG
    printf(
       "   copying N_copy=%d elements (sizeof_N_copy_elements=%d bytes)\n",
           N_copy, (int) sizeof_N_copy_elements);
    #endif
    memcpy(value_buffer, tep->value, sizeof_N_copy_elements);
      }
  }

  return tep->N_elements;
    }
    }
}

/******************************************************************************/

/*
 * This function gets a pointer to a table's header, given the table handle.
 *
 * Arguments:
 * handle = The table handle.
 *
 * Results:
 * If the handle is invalid (i.e. there is no such table), this function
 *    returns NULL.
 * If the handle is valid, this function returns a pointer to the table header.
 */
static
  struct table_header *get_table_header_ptr(int handle)
{
  return ((handle >= 0) && (handle < N_thp_array))
         ? (struct table_header *) thp_array[handle]      /* valid handle */
         : NULL;                                          /* invalid handle */
}

/******************************************************************************/

/*
 * check if key is syntactically "bad" (eg contains '/' character)
 * returns true for bad key, false for ok
 */
static
  bool is_bad_key(const char *key)
{
  assert(key != NULL);

  if (strchr(key, '/') != NULL)
  {
    return true;
  }

  return false;                           /* ok */
}

/******************************************************************************/

/*
 * This function finds the (first) table entry with a given key.
 * Optionally, it also finds the table entry *before* that one in
 * the linked list.
 *
 * Arguments:
 * thp -> The table header.
 * key -> The key to search for.
 * prev_tep_ptr = If this is non-NULL, then this function sets
 *                *prev_tep_ptr to -> the table entry *before* the one
 *                with the given key, or NULL if the table entry with
 *                the given key is the starting one in the table.
 *                Thus if  prev_tep_ptr  is non-NULL, then after this
 *                function returns, the returned result is
 *                (*prev_tep_ptr == NULL) ? thp->head : (*prev_tep)->next
 *
 * Results:
 * The function returns a pointer to the table entry, or NULL if the
 * key isn't found in the table.
 */
static
  struct table_entry *find_table_entry(const struct table_header *thp,
                                       const char *key,
                                       struct table_entry **prev_tep_ptr)
{
  assert(thp != NULL);
  assert(key != NULL);

    {
  const bool case_insensitive_flag
          = thp->flags & UTIL_TABLE_FLAGS_CASE_INSENSITIVE;
  struct table_entry *prev_tep = NULL;
  struct table_entry *tep = thp->head;
  for ( ; tep != NULL ; prev_tep = tep, tep = tep->next)
  {
    if (  case_insensitive_flag
          ?  (Util_StrCmpi(key, tep->key) == 0)
          :  (     strcmp (key, tep->key) == 0)  )
    {
      if (prev_tep_ptr != NULL)
      {
        *prev_tep_ptr = prev_tep;
      }
      return tep;                         /* key found in table */
    }
  }

  return NULL;                            /* key not found in table */
    }
}

/******************************************************************************/

/*@@
  @routine      insert_table_entry
  @desc
                This is an internal function used in implementing
                Util_TableClone() and internal_set().  It allocates
                a new table entry and sets the fields in it to be
                copies of the arguments.
  @enddesc

  @var          thp
  @vtype        struct table_header *
  @vdesc        pointer to the table header
  @endvar

  @var          key
  @vtype        const char *
  @vdesc        pointer to a (C-style null-terminated) string, a copy
                of which is to be the new table entry's key
  @endvar

  @var          type_code
  @vtype        int
  @vdesc        the value to be the new table entry's type code
                (one of the CCTK_VARIABLE_* constants from "cctk_Types.h"),
  @endvar

  @var          N_elements
  @vtype        int (must be >= 0)
  @vdesc        number of elements in  array[]
  @endvar

  @var          array
  @vtype        const T[], where T is one of
                   CCTK_POINTER, CCTK_POINTER_TO_CONST, CCTK_FPOINTER,
                   CCTK_CHAR,
                   CCTK_BYTE,
                   CCTK_INT, CCTK_INT1, CCTK_INT2, CCTK_INT4, CCTK_INT8,
                   CCTK_INT16,
                   CCTK_REAL, CCTK_REAL4, CCTK_REAL8, CCTK_REAL16,
                   CCTK_COMPLEX, CCTK_COMPLEX8, CCTK_COMPLEX16, CCTK_COMPLEX32
                (not all of these may be supported on any given system)
  @vdesc        (pointer to) an array, a copy of which is to be
                the new table entry's value
  @endvar

  @returntype   int
  @returndesc
                0 for ok,<BR>
                UTIL_ERROR_NO_MEMORY            unable to allocate memory,<BR>
                UTIL_ERROR_BAD_INPUT            unknown  type_code
  @endreturndesc
  @@*/
static
  int insert_table_entry(struct table_header *thp,
                         const char *key,
                         int type_code, int N_elements, const void *value)
{
  #ifdef UTIL_TABLE_DEBUG
  printf("insert_table_entry(type_code=%d, N_elements=%d, key=\"%s\")...\n",
         type_code, N_elements, key);
  #endif

    {
  struct table_entry *tep = (struct table_entry *)
                            malloc(sizeof(struct table_entry));
  if (tep == NULL)
  {
    return UTIL_ERROR_NO_MEMORY;            /* can't allocate new table entry */
  }

  tep->key = Util_Strdup(key);
  if (tep->key == NULL)
  {
    free(tep);
    return UTIL_ERROR_NO_MEMORY;         /* can't allocate memory to copy key */
  }

  tep->type_code = type_code;
  tep->N_elements = N_elements;

    {
  const int element_size = CCTK_VarTypeSize(type_code);
  if (element_size < 0)
  {
    free(tep);
    return UTIL_ERROR_BAD_INPUT;        /* unknown  type_code */
  }
  {
  const size_t sizeof_value = N_elements * element_size;
  #ifdef UTIL_TABLE_DEBUG2
  printf("   allocating new buffer of size sizeof_value=%d bytes\n",
         (int) sizeof_value);
  #endif
  void *const buffer = malloc(sizeof_value);
  /*
   * A 0-sized array is (or should be) legal for the table routines.
   * Alas, on some systems malloc(0) returns a NULL pointer, so we must
   * specially check for a 0-sized array in this next test to avoid
   * falsely seeing this as "malloc failed -- we're out of memory".
   */
  if (sizeof_value != 0 && buffer == NULL)
  {
    free(tep->key);
    free(tep);
    return UTIL_ERROR_NO_MEMORY;   /* can't allocate memory for copy of value */
  }
  #ifdef UTIL_TABLE_DEBUG
  printf("   copying sizeof_value=%d bytes into buffer\n", (int) sizeof_value);
  #endif
  memcpy(buffer, value, sizeof_value);
  tep->value = buffer;

  /* insert the table entry into the table's linked list */
  /* (we could insert it anywhere; for simplicity we insert it at the head) */
  tep->next = thp->head;
  thp->head = tep;

  return 0;
    }
    }
    }
}

/******************************************************************************/

/*
 * This function deletes an entry (specified by its key) from a table,
 * freeing its table entry and the pointed-to key and value.
 *
 * Results:
 * The return value is the same as for Util_TableDeleteKey(), i.e.
 *      0 for ok (key existed before this call, and has now been deleted)
 *      -ve for error, including
 *      UTIL_ERROR_TABLE_NO_SUCH_KEY    no such key in table
 */
static
  int delete_table_entry_by_key(struct table_header *thp, const char *key)
{
  struct table_entry *prev_tep;
  struct table_entry *const tep = find_table_entry(thp, key, &prev_tep);
  if (tep == NULL)
  {
    return UTIL_ERROR_TABLE_NO_SUCH_KEY;
  }

  delete_table_entry_by_ptr(thp, prev_tep);
  return 0;                           /* ok: key existed before this call, */
                                      /* and has now been deleted */
}

/******************************************************************************/

/*
 * This function deletes an entry from a table, freeing its table entry
 * and the pointed-to key and value.  The entry to be deleted is specified
 * by a pointer to the *previous* table entry, or NULL to delete the
 * starting entry in the list.
 */
static
  void delete_table_entry_by_ptr(struct table_header *thp,
                                 struct table_entry *prev_tep)
{
/* this is the entry we want to delete */
struct table_entry *const tep = (prev_tep == NULL) ? thp->head
                                                   : prev_tep->next;
assert(tep != NULL);

/* unlink it from the list */
if (prev_tep == NULL)
{
  thp->head = tep->next;
}
else
{
  prev_tep->next = tep->next;
}

assert(tep->key != NULL);
free(tep->key);

assert(tep->N_elements == 0 || tep->value != NULL);
free(tep->value);

free(tep);
}

/******************************************************************************/

/*
 * This function gets a pointer to an iterator, given the iterator handle.
 *
 * Arguments:
 * ihandle = The iterator handle.
 *
 * Results:
 * If the handle is valid, this function returns a pointer to the iterator.
 * If the handle is invalid (i.e. there is no such table), this function
 *    returns NULL.
 */
static
  struct iterator *get_iterator_ptr(int ihandle)
{
  return ((ihandle >= 0) && (ihandle < N_ip_array))
         ? (struct iterator *) ip_array[ihandle]          /* valid handle */
         : NULL;                                          /* invalid handle */
}

/******************************************************************************/

/*
 * This function grows an malloc-allocated array of  void *  pointers
 * via realloc(), initializing the new space to NULL pointers.
 *
 * Arguments:
 * *pN = (in out) array size
 * *pvp_array = (in out) Pointer to growable array of  void *  pointers.
 *
 * Results:
 * This function returns
 *      0 for ok,
 *      -ve for error, including
 *      UTIL_ERROR_NO_MEMORY            can't allocate memory to grow table
 */
static
  int grow_pointer_array(int *pN, void ***pvp_array)
{
const int N = *pN;
  void **vp_array = *pvp_array;
  const int new_N = GROW(N);
  void **new_vp_array = realloc(vp_array, new_N*sizeof(void *));
  if (new_vp_array == NULL)
  {
    return UTIL_ERROR_NO_MEMORY;                          /* can't grow array */
  }

  /* initialize the new space to NULL pointers */
    {
  int i;
  for (i = N ; i < new_N ; ++i)
  {
    new_vp_array[i] = NULL;
  }
    }

  *pvp_array = new_vp_array;
  *pN = new_N;
  return 0;                               /* ok */
}

/******************************************************************************/

/*
 * This function converts the given string into a scalar value of type
 * CCTK_INT or CCTK_REAL.
 *
 * Arguments:
 * string = (in) null-terminated string to be converted into a number
 * scalar = (out) structure defining the type and value of the number
 */
static
  void convert_string_to_number(const char *string, struct scalar_value *scalar)
{
  long int int_scalar;
  double real_scalar;
  char *endptr;


  scalar->datatype = -1;

  if (*string)
  {
    int_scalar = strtol (string, &endptr, 10);
    if (*endptr == 0 && int_scalar != LONG_MIN && int_scalar != LONG_MAX)
    {
      scalar->datatype = CCTK_VARIABLE_INT;
      scalar->value.int_scalar = int_scalar;
    }
    else
    {
      real_scalar = strtod (string, &endptr);
      if (*endptr == 0 && real_scalar != +HUGE_VAL && real_scalar != -HUGE_VAL)
      {
        scalar->datatype = CCTK_VARIABLE_REAL;
        scalar->value.real_scalar = real_scalar;
      }
    }
  }
}

/******************************************************************************/
/***** Table and Iterator Dump Routines ***************************************/
/******************************************************************************/

/*
* This function prints out all the tables and their data structures.
*/
int Util_TablePrintAll(FILE *stream)
{
  int handle;

  fprintf(stream, "N_tables=%d N_thp_array=%d\n", N_tables, N_thp_array);
  for (handle = 0 ; handle < N_thp_array ; ++handle)
  {
    Util_TablePrint(stream, handle);
  }

  return 0;
}

/******************************************************************************/

/*
 * This function prints out a table, giving all internal information.
 */
int Util_TablePrint(FILE *stream, int handle)
{
  fprintf(stream, "thp_array[%d]: ", handle);
    {
  const struct table_header *const thp = get_table_header_ptr(handle);
  if (thp == NULL)
  {
    fprintf(stream, "NULL\n");
  }
  else
  {
    fprintf(stream, "flags=0x%x handle=%d\n", thp->flags, thp->handle);
      {
    const struct table_entry *tep = thp->head;
    for ( ; tep != NULL ; tep = tep->next)
    {
      fprintf(stream, "    [tep=%p]\n", (const void *) tep);
      fprintf(stream, "\tkey=\"%s\"\n", tep->key);
      fprintf(stream, "\ttype_code=%d N_elements=%d\n", tep->type_code, tep->N_elements);
        {
      int i;
      switch  (tep->type_code)
      {
        case CCTK_VARIABLE_BYTE:
          fprintf(stream, "\t[byte]");
            {
          const CCTK_BYTE* const value_ptr_byte = (const CCTK_BYTE*) tep->value;
          for (i = 0 ; i < tep->N_elements ; ++i)
          {
            fprintf(stream, "\t%d", (int) value_ptr_byte[i]);
          }
            }
          break;
        case CCTK_VARIABLE_INT:
          fprintf(stream, "\t[int]");
            {
          const CCTK_INT* const value_ptr_int = (const CCTK_INT*) tep->value;
          for (i = 0 ; i < tep->N_elements ; ++i)
          {
            fprintf(stream, "\t%d", (int) value_ptr_int[i]);
          }
            }
          break;
        case CCTK_VARIABLE_REAL:
          fprintf(stream, "\t[real]");
            {
          const CCTK_REAL* const value_ptr_real = (const CCTK_REAL*) tep->value;
          for (i = 0 ; i < tep->N_elements ; ++i)
          {
            fprintf(stream, "\t%g", (double) value_ptr_real[i]);
          }
            }
          break;
        case CCTK_VARIABLE_CHAR:
          fprintf(stream, "\t[char]");
            {
          const CCTK_CHAR* const value_ptr_char = (const CCTK_CHAR*) tep->value;
          fprintf(stream, "\t\"");
          for (i = 0 ; i < tep->N_elements ; ++i)
          {
            CCTK_CHAR const c = value_ptr_char[i];
            if (c == '"')
              fprintf(stream, "\\\"");
            else if (isprint(c))
              fprintf(stream, "%c", (int) c);
            else
              fprintf(stream, "\\x%02x", (int) c);
          }
          fprintf(stream, "\"");
            }
          break;
        default:
          fprintf(stream, "\t[sorry, don't know how to print this type!]");
          break;
      }
      fprintf(stream, "\n");
        }
    }
      }
  }
    }

  return 0;
}

/******************************************************************************/

/*
 * This function prints out a table into a nice format.
 */
int Util_TablePrintPretty(FILE *stream, int handle)
{
  const struct table_header *const thp = get_table_header_ptr(handle);
  if (! thp)
  {
    fprintf(stream, "NULL\n");
    return 0;
  }
  
  for (const struct table_entry *tep = thp->head; tep; tep = tep->next)
  {
    if (tep != thp->head)
    {
      fprintf(stream, " ");
    }
    fprintf(stream, "%s=", tep->key);
    switch  (tep->type_code)
    {
    case CCTK_VARIABLE_BYTE:
      {
        const CCTK_BYTE* const value_ptr_byte = (const CCTK_BYTE*)tep->value;
        if (tep->N_elements != 1)
        {
          fprintf(stream, "{");
        }
        for (int i = 0; i < tep->N_elements; ++i)
        {
          if (i != 0)
          {
            fprintf(stream, " ");
          }
          fprintf(stream, "0x%02x", (int)value_ptr_byte[i]);
        }
        if (tep->N_elements != 1)
        {
          fprintf(stream, "}");
        }
      }
      break;
    case CCTK_VARIABLE_INT:
      {
        const CCTK_INT* const value_ptr_int = (const CCTK_INT*)tep->value;
        if (tep->N_elements != 1)
        {
          fprintf(stream, "{");
        }
        for (int i = 0; i < tep->N_elements; ++i)
        {
          if (i != 0)
          {
            fprintf(stream, " ");
          }
          fprintf(stream, "%d", (int)value_ptr_int[i]);
        }
        if (tep->N_elements != 1)
        {
          fprintf(stream, "}");
        }
      }
      break;
    case CCTK_VARIABLE_REAL:
      {
        const CCTK_REAL* const value_ptr_real = (const CCTK_REAL*)tep->value;
        if (tep->N_elements != 1)
        {
          fprintf(stream, "{");
        }
        for (int i = 0; i < tep->N_elements; ++i)
        {
          if (i != 0)
          {
            fprintf(stream, " ");
          }
          fprintf(stream, "%#.17g", (double)value_ptr_real[i]);
        }
        if (tep->N_elements != 1)
        {
          fprintf(stream, "}");
        }
      }
      break;
    case CCTK_VARIABLE_CHAR:
      {
        const CCTK_CHAR* const value_ptr_char = (const CCTK_CHAR*)tep->value;
        fprintf(stream, "\"");
        for (int i = 0; i < tep->N_elements; ++i)
        {
          CCTK_CHAR const c = value_ptr_char[i];
          if (c == '"')
            fprintf(stream, "\\\"");
          else if (isprint(c))
            fprintf(stream, "%c", (int)c);
          else
            fprintf(stream, "\\x%02x", (int)(unsigned char)c);
        }
        fprintf(stream, "\"");
      }
      break;
    default:
      fprintf(stream, "[unsupported type]");
    }
  }
  return 0;
}

/******************************************************************************/

/*
 * This function prints out all the iterators and their data structures.
 */
int Util_TablePrintAllIterators(FILE *stream)
{
  int ihandle;

  fprintf(stream, "N_iterators=%d N_ip_array=%d\n", N_iterators, N_ip_array);
  for (ihandle = 0 ; ihandle < N_ip_array ; ++ihandle)
  {
    const struct iterator *const ip = get_iterator_ptr(ihandle);
    fprintf(stream, "ip_array[%d]: ", ihandle);
    if (ip == NULL)
    {
      fprintf(stream, "NULL\n");
    }
    else
    {
      fprintf(stream, "thp=%p tep=%p\n", (const void *) ip->thp, (const void *) ip->tep);
    }
  }

  return 0;
}
