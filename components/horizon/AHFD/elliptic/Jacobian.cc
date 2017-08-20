// Jacobian.cc -- generic routines for the Jacobian matrix
// $Header$

//
// <<<literal contents of "lapack.hh">>>
//
// decode_Jacobian_store_solve_method -- decode string into internal enum
// new_Jacobian -- object factory for Jacobian objects
//

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "../jtutil/util_Table.h"
#include "../AHFD_macros.h"


#include "../jtutil/util.hh"
#include "../jtutil/array.hh"
#include "../jtutil/cpm_map.hh"
#include "../jtutil/linear_map.hh"

#include "../patch/coords.hh"
#include "../patch/grid.hh"
#include "../patch/fd_grid.hh"
#include "../patch/patch.hh"
#include "../patch/patch_edge.hh"
#include "../patch/patch_interp.hh"
#include "../patch/ghost_zone.hh"
#include "../patch/patch_system.hh"

#include "Jacobian.hh"
#include "dense_Jacobian.hh"
#include "row_sparse_Jacobian.hh"

// all the code in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function decodes a character string specifying a specific type
// (derived class) of Jacobian, into an internal enum.
//
enum Jacobian_store_solve_method
  decode_Jacobian_store_solve_method
	(const char Jacobian_store_solve_method_string[])
{
if	(STRING_EQUAL(Jacobian_store_solve_method_string,
		      "dense matrix/LAPACK"))
   then {
  #ifdef HAVE_DENSE_JACOBIAN__LAPACK
	return Jacobian__dense_matrix__LAPACK;
  #endif
	}

else if (STRING_EQUAL(Jacobian_store_solve_method_string,
		      "row-oriented sparse matrix/ILUCG"))
   then {
  #ifdef HAVE_ROW_SPARSE_JACOBIAN__ILUCG
	return Jacobian__row_sparse_matrix__ILUCG;
  #endif
	}

else if (STRING_EQUAL(Jacobian_store_solve_method_string,
		      "row-oriented sparse matrix/UMFPACK"))
   then {
  #ifdef HAVE_ROW_SPARSE_JACOBIAN__UMFPACK
	return Jacobian__row_sparse_matrix__UMFPACK;
  #endif
	}

else	error_exit(ERROR_EXIT,
"decode_Jacobian_store_solve_method():\n"
"        unknown Jacobian_store_solve_method_string=\"%s\"!\n",
		   Jacobian_store_solve_method_string);		/*NOTREACHED*/

// fall through to here ==> we recognize the matrix store_solve_method,
//                          but it's not configured
error_exit(ERROR_EXIT,
"\n"
"   decode_Jacobian_store_solve_method():\n"
"        Jacobian store_solve_method=\"%s\"\n"
"        is not configured in this binary; see \"src/include/config.hh\"\n"
"        for details on what methods are configured and how to change this\n"
	   ,
	   Jacobian_store_solve_method_string);			/*NOTREACHED*/
}

//******************************************************************************

//
// This function is an "object factory" for Jacobian objects: it constructs
// and returns a pointer to a new-allocated Jacobian object of the
// specified derived type.
//
// FIXME: the patch system shouldn't really have to be non-const, but
//	  the Jacobian constructors all require this to allow the
//	  linear solvers to directly update gridfns
//
Jacobian* new_Jacobian(enum Jacobian_store_solve_method Jac_method,
		       patch_system& ps,
		       bool print_msg_flag /* = false */)
{
switch	(Jac_method)
	{
#ifdef HAVE_DENSE_JACOBIAN__LAPACK
  case Jacobian__dense_matrix__LAPACK:
	return new dense_Jacobian__LAPACK(ps, print_msg_flag);
#endif

#ifdef HAVE_ROW_SPARSE_JACOBIAN__ILUCG
  case Jacobian__row_sparse_matrix__ILUCG:
	return new row_sparse_Jacobian__ILUCG(ps, print_msg_flag);
#endif

#ifdef HAVE_ROW_SPARSE_JACOBIAN__UMFPACK
  case Jacobian__row_sparse_matrix__UMFPACK:
	return new row_sparse_Jacobian__UMFPACK(ps, print_msg_flag);
#endif

  default:
	error_exit(ERROR_EXIT,
		   "new_Jacobian(): unknown method=(int)%d!\n",
		   int(Jac_method));				/*NOTREACHED*/
	}
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

	  }	// namespace AHFinderDirect
