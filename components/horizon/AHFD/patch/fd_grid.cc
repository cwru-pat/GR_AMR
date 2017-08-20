// fd_grid.cc -- grid with finite differencing operations
// $Header$
//
// fd_grid::dx_coeff
// fd_grid::dxx_coeff
//

#include <stdio.h>
#include <assert.h>
#include <math.h>


  #include "../cctk.h"

#include "../jtutil/util.hh"
#include "../jtutil/array.hh"
#include "../jtutil/linear_map.hh"

#include "coords.hh"
#include "grid.hh"
#include "fd_grid.hh"

// all the code in this file is inside this namespace
namespace AHFD
	  {

//*****************************************************************************

//
// This function computes a single coefficient of a 1st derivative
// molecule, for unit grid spacing.
//
//static
  fp fd_grid::dx_coeff(int m)
{
#ifndef FINITE_DIFF_ORDER
  #error "must define FINITE_DIFF_ORDER!"
#endif

switch	(m)
	{
#if   (FINITE_DIFF_ORDER == 2)
  case -1:	return FD_GRID__ORDER2__DX__COEFF_M1;
  case  0:	return FD_GRID__ORDER2__DX__COEFF_0;
  case +1:	return FD_GRID__ORDER2__DX__COEFF_P1;
#elif (FINITE_DIFF_ORDER == 4)
  case -2:	return FD_GRID__ORDER4__DX__COEFF_M2;
  case -1:	return FD_GRID__ORDER4__DX__COEFF_M1;
  case  0:	return FD_GRID__ORDER4__DX__COEFF_0;
  case +1:	return FD_GRID__ORDER4__DX__COEFF_P1;
  case +2:	return FD_GRID__ORDER4__DX__COEFF_P2;
#else
  #error "unknown value " FINITE_DIFF_ORDER " for FINITE_DIFF_ORDER!"
#endif
default:
	error_exit(ERROR_EXIT,
"***** fd_grid::dx_coeff(): m=%d is outside order=%d molecule radius=%d\n",
		   m, FINITE_DIFF_ORDER, FD_GRID__MOL_RADIUS);	/*NOTREACHED*/
	}
}

//*****************************************************************************

//
// This function computes a single coefficient of a 2nd derivative
// molecule, for unit grid spacing.
//
//static
  fp fd_grid::dxx_coeff(int m)
{
#ifndef FINITE_DIFF_ORDER
  #error "must define FINITE_DIFF_ORDER!"
#endif

switch	(m)
	{
#if   (FINITE_DIFF_ORDER == 2)
  case -1:	return FD_GRID__ORDER2__DXX__COEFF_M1;
  case  0:	return FD_GRID__ORDER2__DXX__COEFF_0;
  case +1:	return FD_GRID__ORDER2__DXX__COEFF_P1;
#elif (FINITE_DIFF_ORDER == 4)
  case -2:	return FD_GRID__ORDER4__DXX__COEFF_M2;
  case -1:	return FD_GRID__ORDER4__DXX__COEFF_M1;
  case  0:	return FD_GRID__ORDER4__DXX__COEFF_0;
  case +1:	return FD_GRID__ORDER4__DXX__COEFF_P1;
  case +2:	return FD_GRID__ORDER4__DXX__COEFF_P2;
#else
  #error "unknown value " FINITE_DIFF_ORDER " for FINITE_DIFF_ORDER!"
#endif
default:
	error_exit(ERROR_EXIT,
"***** fd_grid::dxx_coeff(): m=%d is outside order=%d molecule radius=%d\n",
		   m, FINITE_DIFF_ORDER, FD_GRID__MOL_RADIUS);	/*NOTREACHED*/
	}
}

//******************************************************************************

	  }	// namespace AHFinderDirect
