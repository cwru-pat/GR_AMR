#ifndef COSMO_UTILS_MATH_H
#define COSMO_UTILS_MATH_H

#include "../cosmo_types.h"
//#include "../cosmo_includes.h"
#include "SAMRAI/pdat/MDA_Access.h"

namespace cosmo
{

/**
 * @brief Kreiss-Oliger dissipation for 2nd-order stencils
 * @details
 * Works only for 2nd-order accurate (Odx2) derivatives
 * accuracy a = 2r - 2
 * r = (a + 2)/2
 * for a = 2, r = 2
 * Q = (-1)^r * (dx)^(2r-1) D_+^r D_-^r / 2^(2r)
 *   = dx^3 / 2^6 D_+^2 D_-^2
 *   = dx^3 / 64 * [stencil: 1, -4, 6, -4, 1]
 * 
 * for a = 8, r = 5
 * Q = (-1)^r * (dx)^(2r-1) D_+^r D_-^r / 2^(2r)
 *   = -dx^9 / 1024 * D_+^5 D_-^5
 *   = -dx^9 / 1024 * [stencil: 1, -10, 45, -120, 210, -252, 210, -120, 45, -10, 1]
 *   
 * TODO: other orders?
 * 
 * @param i gridpoint in x-dir
 * @param j gridpoint in y-dir
 * @param k gridpoint in z-dir
 * @param field releavnt field
 * @return dissipation factor
 */
inline real_t KO_dissipation_Q(
  idx_t i, idx_t j, idx_t k, arr_t & field, const double dx[], real_t ko_coeff)
{
  if(ko_coeff == 0)
    return 0;

# if STENCIL_ORDER == 2
    real_t stencil = (
        1.0*field(i-2,j,k) + 1.0*field(i,j-2,k) + 1.0*field(i,j,k-2)
      - 4.0*field(i-1,j,k) - 4.0*field(i,j-1,k) - 4.0*field(i,j,k-1)
      + 6.0*field(i  ,j,k) + 6.0*field(i,j  ,k) + 6.0*field(i,j,k  )
      - 4.0*field(i+1,j,k) - 4.0*field(i,j+1,k) - 4.0*field(i,j,k+1)
      + 1.0*field(i+2,j,k) + 1.0*field(i,j+2,k) + 1.0*field(i,j,k+2)
    )/pow(dx, 4.0);
    real_t dissipation = ko_coeff*pow(dx, 3.0)/64.0*stencil;
    return dissipation;
# endif

# if STENCIL_ORDER == 4
    real_t stencil = (
      1.0*field(i+3,j,k) 
      - 6.0*field(i+2,j,k)
      + 15.0*field(i+1,j,k) 
      - 20.0*field(i  ,j,k) 
      + 15.0*field(i-1,j,k) 
      - 6.0*field(i-2,j,k) 
      + 1.0*field(i-3,j,k) 
    )/dx[0]
      + (1.0*field(i,j+3,k) 
      - 6.0*field(i,j+2,k)
      + 15.0*field(i,j+1,k) 
      - 20.0*field(i  ,j,k) 
      + 15.0*field(i,j-1,k) 
      - 6.0*field(i,j-2,k) 
      + 1.0*field(i,j-3,k) 
    )/dx[1]
      + (1.0*field(i,j,k+3) 
      - 6.0*field(i,j,k+2)
      + 15.0*field(i,j,k+1) 
      - 20.0*field(i  ,j,k) 
      + 15.0*field(i,j,k-1) 
      - 6.0*field(i,j,k-2) 
      + 1.0*field(i,j,k-3) 
      )/dx[2];

    return stencil * ko_coeff /64.0;
    
# endif

    
# if STENCIL_ORDER == 8
    real_t stencil = (
          1.0*field(i-5,j,k) +   1.0*field(i,j-5,k) +   1.0*field(i,j,k-5)
      -  10.0*field(i-4,j,k) -  10.0*field(i,j-4,k) -  10.0*field(i,j,k-4)
      +  45.0*field(i-3,j,k) +  45.0*field(i,j-3,k) +  45.0*field(i,j,k-3)
      - 120.0*field(i-2,j,k) - 120.0*field(i,j-2,k) - 120.0*field(i,j,k-2)
      + 210.0*field(i-1,j,k) + 210.0*field(i,j-1,k) + 210.0*field(i,j,k-1)
      - 252.0*field(i  ,j,k) - 252.0*field(i,j  ,k) - 252.0*field(i,j,k  )
      + 210.0*field(i+1,j,k) + 210.0*field(i,j+1,k) + 210.0*field(i,j,k+1)
      - 120.0*field(i+2,j,k) - 120.0*field(i,j+2,k) - 120.0*field(i,j,k+2)
      +  45.0*field(i+3,j,k) +  45.0*field(i,j+3,k) +  45.0*field(i,j,k+3)
      -  10.0*field(i+4,j,k) -  10.0*field(i,j+4,k) -  10.0*field(i,j,k+4)
      +   1.0*field(i+5,j,k) +   1.0*field(i,j+5,k) +   1.0*field(i,j,k+5)
    )/pow(dx, 10.0);
    real_t dissipation = -ko_coeff*pow(dx, 9.0)/1024.0*stencil;
    return dissipation;
# endif

  return 0.0;
}

inline real_t derivative_Odx2(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return ((
        - 1.0/2.0*field(i-1,j,k)
        + 1.0/2.0*field(i+1,j,k)
      )/dx[0]);
      break;
    case 2:
      return ((
        - 1.0/2.0*field(i,j-1,k)
        + 1.0/2.0*field(i,j+1,k)
      )/dx[1]);
      break;
    case 3:
      return ((
        - 1.0/2.0*field(i,j,k-1)
        + 1.0/2.0*field(i,j,k+1)
      )/dx[2]);
      break;
  }

  /* XXX */
  return 0;
}

inline real_t forward_derivative_Odx2(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return ((
        - 1.0/2.0*field(i+2,j,k)
        + 2.0*field(i+1,j,k)
        - 3.0/2.0*field(i,j,k)
      )/dx[0]);
      break;
    case 2:
      return ((
        - 1.0/2.0*field(i,j+2,k)
        + 2.0*field(i,j+1,k)
        - 3.0/2.0*field(i,j,k)
      )/dx[1]);
      break;
    case 3:
      return ((
        - 1.0/2.0*field(i,j,k+2)
        + 2.0*field(i,j,k+1)
        - 3.0/2.0*field(i,j,k)
      )/dx[2]);
      break;
  }

  /* XXX */
  return 0;
}



inline real_t backward_derivative_Odx2(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return ((
        + 1.0/2.0*field(i-2,j,k)
        - 2.0*field(i-1,j,k)
        + 3.0/2.0*field(i,j,k)
      )/dx[0]);
      break;
    case 2:
      return ((
        + 1.0/2.0*field(i,j-2,k)
        - 2.0*field(i,j-1,k)
        + 3.0/2.0*field(i,j,k)
      )/dx[1]);
      break;
    case 3:
      return ((
        + 1.0/2.0*field(i,j,k-2)
        - 2.0*field(i,j,k-1)
        + 3.0/2.0*field(i,j,k)
      )/dx[2]);
      break;
  }

  /* XXX */
  return 0;
}

 
inline real_t derivative_Odx4(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
        + 1.0/12.0*field(i-2,j,k)
        - 2.0/3.0*field(i-1,j,k)
        + 2.0/3.0*field(i+1,j,k)
        - 1.0/12.0*field(i+2,j,k)
      )/dx[0];
      break;
    case 2:
      return (
        + 1.0/12.0*field(i,j-2,k)
        - 2.0/3.0*field(i,j-1,k)
        + 2.0/3.0*field(i,j+1,k)
        - 1.0/12.0*field(i,j+2,k)
      )/dx[1];
      break;
    case 3:
      return (
        + 1.0/12.0*field(i,j,k-2)
        - 2.0/3.0*field(i,j,k-1)
        + 2.0/3.0*field(i,j,k+1)
        - 1.0/12.0*field(i,j,k+2)
      )/dx[2];
      break;
  }

  /* XXX */
  return 0;
}

inline real_t forward_derivative_Odx4(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
        - 1.0/4.0*field(i+4,j,k)
        + 4.0/3.0*field(i+3,j,k)
        - 3.0*field(i+2,j,k)
        + 4.0*field(i+1,j,k)
        - 25.0/12.0*field(i,j,k)
      )/dx[0];
      break;
    case 2:
      return (
        - 1.0/4.0*field(i,j+4,k)
        + 4.0/3.0*field(i,j+3,k)
        - 3.0*field(i,j+2,k)
        + 4.0*field(i,j+1,k)
        - 25.0/12.0*field(i,j,k)
      )/dx[1];
      break;
    case 3:
      return (
        - 1.0/4.0*field(i,j,k+4)
        + 4.0/3.0*field(i,j,k+3)
        - 3.0*field(i,j,k+2)
        + 4.0*field(i,j,k+1)
        - 25.0/12.0*field(i,j,k)
      )/dx[2];
      break;
  }

  /* XXX */
  return 0;
}

inline real_t backward_derivative_Odx4(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
        + 1.0/4.0*field(i-4,j,k)
        - 4.0/3.0*field(i-3,j,k)
        + 3.0*field(i-2,j,k)
        - 4.0*field(i-1,j,k)
        + 25.0/12.0*field(i,j,k)
      )/dx[0];
      break;
    case 2:
      return (
        + 1.0/4.0*field(i,j-4,k)
        - 4.0/3.0*field(i,j-3,k)
        + 3.0*field(i,j-2,k)
        - 4.0*field(i,j-1,k)
        + 25.0/12.0*field(i,j,k)
      )/dx[1];
      break;
    case 3:
      return (
        + 1.0/4.0*field(i,j,k-4)
        - 4.0/3.0*field(i,j,k-3)
        + 3.0*field(i,j,k-2)
        - 4.0*field(i,j,k-1)
        + 25.0/12.0*field(i,j,k)
      )/dx[2];
      break;
  }
  return 0;
}
 
inline real_t forward_derivative_Odx6(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
        - 49.0/20.0*field(i,j,k)
        + 6.0*field(i+1,j,k)
        - 15.0/2.0*field(i+2,j,k)
        + 20.0/3.0*field(i+3,j,k)
        - 15.0/4.0*field(i+4,j,k)
        + 6.0/5.0*field(i+5,j,k)
        - 1.0/6.0*field(i+6,j,k)
      )/dx[0];
      break;
    case 2:
      return (
        - 49.0/20.0*field(i,j,k)
        + 6.0*field(i,j+1,k)
        - 15.0/2.0*field(i,j+2,k)
        + 20.0/3.0*field(i,j+3,k)
        - 15.0/4.0*field(i,j+4,k)
        + 6.0/5.0*field(i,j+5,k)
        - 1.0/6.0*field(i,j+6,k)
      )/dx[1];
      break;
    case 3:
      return (
        - 49.0/20.0*field(i,j,k)
        + 6.0*field(i,j,k+1)
        - 15.0/2.0*field(i,j,k+2)
        + 20.0/3.0*field(i,j,k+3)
        - 15.0/4.0*field(i,j,k+4)
        + 6.0/5.0*field(i,j,k+5)
        - 1.0/6.0*field(i,j,k+6)
      )/dx[2];
      break;
  }

  /* XXX */
  return 0;
}

inline real_t backward_derivative_Odx6(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
        + 49.0/20.0*field(i,j,k)
        - 6.0*field(i-1,j,k)
        + 15.0/2.0*field(i-2,j,k)
        - 20.0/3.0*field(i-3,j,k)
        + 15.0/4.0*field(i-4,j,k)
        - 6.0/5.0*field(i-5,j,k)
        + 1.0/6.0*field(i-6,j,k)
      )/dx[0];
      break;
    case 2:
      return (
        + 49.0/20.0*field(i,j,k)
        - 6.0*field(i,j-1,k)
        + 15.0/2.0*field(i,j-2,k)
        - 20.0/3.0*field(i,j-3,k)
        + 15.0/4.0*field(i,j-4,k)
        - 6.0/5.0*field(i,j-5,k)
        + 1.0/6.0*field(i,j-6,k)
      )/dx[1];
      break;
    case 3:
      return (
        + 49.0/20.0*field(i,j,k)
        - 6.0*field(i,j,k-1)
        + 15.0/2.0*field(i,j,k-2)
        - 20.0/3.0*field(i,j,k-3)
        + 15.0/4.0*field(i,j,k-4)
        - 6.0/5.0*field(i,j,k-5)
        + 1.0/6.0*field(i,j,k-6)
      )/dx[2];
      break;
  }

  /* XXX */
  return 0;
}

inline real_t derivative_Odx6(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
        - 1.0/60.0*field(i-3,j,k)
        + 3.0/20.0*field(i-2,j,k)
        - 3.0/4.0*field(i-1,j,k)
        + 3.0/4.0*field(i+1,j,k)
        - 3.0/20.0*field(i+2,j,k)
        + 1.0/60.0*field(i+3,j,k)
      )/dx[0];
      break;
    case 2:
      return (
        - 1.0/60.0*field(i,j-3,k)
        + 3.0/20.0*field(i,j-2,k)
        - 3.0/4.0*field(i,j-1,k)
        + 3.0/4.0*field(i,j+1,k)
        - 3.0/20.0*field(i,j+2,k)
        + 1.0/60.0*field(i,j+3,k)
      )/dx[1];
      break;
    case 3:
      return (
        - 1.0/60.0*field(i,j,k-3)
        + 3.0/20.0*field(i,j,k-2)
        - 3.0/4.0*field(i,j,k-1)
        + 3.0/4.0*field(i,j,k+1)
        - 3.0/20.0*field(i,j,k+2)
        + 1.0/60.0*field(i,j,k+3)
      )/dx[2];
      break;
  }

  /* XXX */
  return 0;
}

 
inline real_t derivative_Odx8(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
        ( 1.0/280.0*field(i-4,j,k) - 1.0/280.0*field(i+4,j,k) )
        - ( 4.0/105.0*field(i-3,j,k) - 4.0/105.0*field(i+3,j,k) )
        + ( 1.0/5.0*field(i-2,j,k) - 1.0/5.0*field(i+2,j,k) )
        - ( 4.0/5.0*field(i-1,j,k) - 4.0/5.0*field(i+1,j,k) )
      )/dx[0];
      break;
    case 2:
      return (
        ( 1.0/280.0*field(i,j-4,k) - 1.0/280.0*field(i,j+4,k) )
        - ( 4.0/105.0*field(i,j-3,k) - 4.0/105.0*field(i,j+3,k) )
        + ( 1.0/5.0*field(i,j-2,k) - 1.0/5.0*field(i,j+2,k) )
        - ( 4.0/5.0*field(i,j-1,k) - 4.0/5.0*field(i,j+1,k) )
      )/dx[1];
      break;
    case 3:
      return (
        ( 1.0/280.0*field(i,j,k-4) - 1.0/280.0*field(i,j,k+4) )
        - ( 4.0/105.0*field(i,j,k-3) - 4.0/105.0*field(i,j,k+3) )
        + ( 1.0/5.0*field(i,j,k-2) - 1.0/5.0*field(i,j,k+2) )
        - ( 4.0/5.0*field(i,j,k-1) - 4.0/5.0*field(i,j,k+1) )
      )/dx[2];
      break;
  }

  /* XXX */
  return 0;
}

inline real_t forward_derivative_Odx8(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
        - 761.0/280.0*field(i,j,k) + 8.0*field(i+1,j,k) 
        - 14.0*field(i+2,j,k) + 56.0/3.0*field(i+3,j,k) 
        - 35.0/2.0*field(i+4,j,k) + 56.0/5.0*field(i+5,j,k) 
        - 14.0/3.0*field(i+6,j,k) + 8.0/7.0*field(i+7,j,k)
        - 1.0/8.0*field(i+8,j,k)
      )/dx[0];
      break;
    case 2:
      return (
        - 761.0/280.0*field(i,j,k) + 8.0*field(i,j+1,k) 
        - 14.0*field(i,j+2,k) + 56.0/3.0*field(i,j+3,k) 
        - 35.0/2.0*field(i,j+4,k) + 56.0/5.0*field(i,j+5,k) 
        - 14.0/3.0*field(i,j+6,k) + 8.0/7.0*field(i,j+7,k)
        - 1.0/8.0*field(i,j+8,k)
      )/dx[1];
      break;
    case 3:
      return (
        - 761.0/280.0*field(i,j,k) + 8.0*field(i,j,k+1) 
        - 14.0*field(i,j,k+2) + 56.0/3.0*field(i,j,k+3) 
        - 35.0/2.0*field(i,j,k+4) + 56.0/5.0*field(i,j,k+5) 
        - 14.0/3.0*field(i,j,k+6) + 8.0/7.0*field(i,j,k+7)
        - 1.0/8.0*field(i,j,k+8)
      )/dx[2];
      break;
  }

  /* XXX */
  return 0;
}

inline real_t backward_derivative_Odx8(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
        + 761.0/280.0*field(i,j,k) - 8.0*field(i-1,j,k) 
        + 14.0*field(i-2,j,k) - 56.0/3.0*field(i-3,j,k) 
        + 35.0/2.0*field(i-4,j,k) - 56.0/5.0*field(i-5,j,k) 
        + 14.0/3.0*field(i-6,j,k) - 8.0/7.0*field(i-7,j,k)
        + 1.0/8.0*field(i-8,j,k)
      )/dx[0];
      break;
    case 2:
      return (
        + 761.0/280.0*field(i,j,k) - 8.0*field(i,j-1,k) 
        + 14.0*field(i,j-2,k) - 56.0/3.0*field(i,j-3,k) 
        + 35.0/2.0*field(i,j-4,k) - 56.0/5.0*field(i,j-5,k) 
        + 14.0/3.0*field(i,j-6,k) - 8.0/7.0*field(i,j-7,k)
        + 1.0/8.0*field(i,j-8,k)
      )/dx[1];
      break;
    case 3:
      return (
        + 761.0/280.0*field(i,j,k) - 8.0*field(i,j,k-1) 
        + 14.0*field(i,j,k-2) - 56.0/3.0*field(i,j,k-3) 
        + 35.0/2.0*field(i,j,k-4) - 56.0/5.0*field(i,j,k-5) 
        + 14.0/3.0*field(i,j,k-6) - 8.0/7.0*field(i,j,k-7)
        + 1.0/8.0*field(i,j,k-8)
      )/dx[2];
      break;
  }

  /* XXX */
  return 0;
}

 
inline real_t mixed_derivative_stencil_Odx2(
  idx_t i, idx_t j, idx_t k, int d1, int d2, arr_t & field, const double dx[])
{
  if( (d1 == 1 && d2 == 2) || (d1 == 2 && d2 == 1) ) {
    return (
      - field(i+1,j-1,k) + field(i+1,j+1,k)
      + field(i-1,j-1,k) - field(i-1,j+1,k)
    )/4.0/dx[0]/dx[1];
  }

  if( (d1 == 1 && d2 == 3) || (d1 == 3 && d2 == 1) ) {
    return (
      - field(i+1,j,k-1) + field(i+1,j,k+1)
      + field(i-1,j,k-1) - field(i-1,j,k+1)
    )/4.0/dx[0]/dx[2];
  }

  if( (d1 == 3 && d2 == 2) || (d1 == 2 && d2 == 3) ) {
    return (
      - field(i,j+1,k-1) + field(i,j+1,k+1)
      + field(i,j-1,k-1) - field(i,j-1,k+1)
    )/4.0/dx[1]/dx[2];
  }

  /* XXX */
  return 0;
}

inline real_t mixed_derivative_stencil_Odx4(
  idx_t i, idx_t j, idx_t k, int d1, int d2, arr_t & field, const double dx[])
{
  if( (d1 == 1 && d2 == 2) || (d1 == 2 && d2 == 1) ) {
    return (
      (
        - field(i+1,j-1,k) + field(i+1,j+1,k)
        + field(i-1,j-1,k) - field(i-1,j+1,k)
      ) - 1.0/16.0*(
        - field(i+2,j-2,k) + field(i+2,j+2,k)
        + field(i-2,j-2,k) - field(i-2,j+2,k)
      )
    )/3.0/dx[0]/dx[1];
  }

  if( (d1 == 1 && d2 == 3) || (d1 == 3 && d2 == 1) ) {
    return (
      (
        - field(i+1,j,k-1) + field(i+1,j,k+1)
        + field(i-1,j,k-1) - field(i-1,j,k+1)
      ) - 1.0/16.0*(
        - field(i+2,j,k-2) + field(i+2,j,k+2)
        + field(i-2,j,k-2) - field(i-2,j,k+2)
      )
    )/3.0/dx[0]/dx[2];
  }

  if( (d1 == 3 && d2 == 2) || (d1 == 2 && d2 == 3) ) {
    return (
      (
        - field(i,j+1,k-1) + field(i,j+1,k+1)
        + field(i,j-1,k-1) - field(i,j-1,k+1)
      ) - 1.0/16.0*(
        - field(i,j+2,k-2) + field(i,j+2,k+2)
        + field(i,j-2,k-2) - field(i,j-2,k+2)
      )
    )/3.0/dx[1]/dx[2];
  }

  /* XXX */
  return 0;
}

inline real_t mixed_derivative_stencil_Odx6(idx_t i, idx_t j, idx_t k, int d1, int d2, arr_t & field, const double dx[])
{
  if( (d1 == 1 && d2 == 2) || (d1 == 2 && d2 == 1) ) {
    return (
      135.0*(
        - field(i+1,j-1,k) + field(i+1,j+1,k)
        + field(i-1,j-1,k) - field(i-1,j+1,k)
      ) - 27.0/2.0*(
        - field(i+2,j-2,k) + field(i+2,j+2,k)
        + field(i-2,j-2,k) - field(i-2,j+2,k)
      ) + (
        - field(i+3,j-3,k) + field(i+3,j+3,k)
        + field(i-3,j-3,k) - field(i-3,j+3,k)
      )
    )/360.0/dx[0]/dx[1];
  }

  if( (d1 == 1 && d2 == 3) || (d1 == 3 && d2 == 1) ) {
    return (
      135.0*(
        - field(i+1,j,k-1) + field(i+1,j,k+1)
        + field(i-1,j,k-1) - field(i-1,j,k+1)
      ) - 27.0/2.0*(
        - field(i+2,j,k-2) + field(i+2,j,k+2)
        + field(i-2,j,k-2) - field(i-2,j,k+2)
      ) + (
        - field(i+3,j,k-3) + field(i+3,j,k+3)
        + field(i-3,j,k-3) - field(i-3,j,k+3)
      )
    )/360.0/dx[0]/dx[2];
  }

  if( (d1 == 3 && d2 == 2) || (d1 == 2 && d2 == 3) ) {
    return (
      135.0*(
        - field(i,j+1,k-1) + field(i,j+1,k+1)
        + field(i,j-1,k-1) - field(i,j-1,k+1)
      ) - 27.0/2.0*(
        - field(i,j+2,k-2) + field(i,j+2,k+2)
        + field(i,j-2,k-2) - field(i,j-2,k+2)
      ) + (
        - field(i,j+3,k-3) + field(i,j+3,k+3)
        + field(i,j-3,k-3) - field(i,j-3,k+3)
      )
    )/360.0/dx[1]/dx[2];
  }

  /* XXX */
  return 0;
}

inline real_t mixed_derivative_stencil_Odx8(idx_t i, idx_t j, idx_t k, int d1,
 int d2, arr_t & field, const double dx[])
{
  if( (d1 == 1 && d2 == 2) || (d1 == 2 && d2 == 1) ) {
    return (
      2.0/5.0*(
        - field(i+1,j-1,k) + field(i+1,j+1,k)
        + field(i-1,j-1,k) - field(i-1,j+1,k)
      ) - 1.0/20.0*(
        - field(i+2,j-2,k) + field(i+2,j+2,k)
        + field(i-2,j-2,k) - field(i-2,j+2,k)
      ) + 2.0/315.0*(
        - field(i+3,j-3,k) + field(i+3,j+3,k)
        + field(i-3,j-3,k) - field(i-3,j+3,k)
      ) - 1.0/2240.0*(
        - field(i+4,j-4,k) + field(i+4,j+4,k)
        + field(i-4,j-4,k) - field(i-4,j+4,k)
      )
    )/dx[0]/dx[1];
  }

  if( (d1 == 1 && d2 == 3) || (d1 == 3 && d2 == 1) ) {
    return (
      2.0/5.0*(
        - field(i+1,j,k-1) + field(i+1,j,k+1)
        + field(i-1,j,k-1) - field(i-1,j,k+1)
      ) - 1.0/20.0*(
        - field(i+2,j,k-2) + field(i+2,j,k+2)
        + field(i-2,j,k-2) - field(i-2,j,k+2)
      ) + 2.0/315.0*(
        - field(i+3,j,k-3) + field(i+3,j,k+3)
        + field(i-3,j,k-3) - field(i-3,j,k+3)
      ) - 1.0/2240.0*(
        - field(i+4,j,k-4) + field(i+4,j,k+4)
        + field(i-4,j,k-4) - field(i-4,j,k+4)
      )
    )/dx[0]/dx[2];
  }

  if( (d1 == 3 && d2 == 2) || (d1 == 2 && d2 == 3) ) {
    return (
      2.0/5.0*(
        - field(i,j+1,k-1) + field(i,j+1,k+1)
        + field(i,j-1,k-1) - field(i,j-1,k+1)
      ) - 1.0/20.0*(
        - field(i,j+2,k-2) + field(i,j+2,k+2)
        + field(i,j-2,k-2) - field(i,j-2,k+2)
      ) + 2.0/315.0*(
        - field(i,j+3,k-3) + field(i,j+3,k+3)
        + field(i,j-3,k-3) - field(i,j-3,k+3)
      ) - 1.0/2240.0*(
        - field(i,j+4,k-4) + field(i,j+4,k+4)
        + field(i,j-4,k-4) - field(i,j-4,k+4)
      )
    )/dx[1]/dx[2];
  }

  /* XXX */
  return 0;
}

inline real_t double_derivative_stencil_Odx2(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
          field(i-1,j,k)
          - 2.0*field(i-0,j,k)
          + field(i+1,j,k)
        )/dx[0]/dx[0];
      break;
    case 2:
      return (
          field(i,j-1,k)
          - 2.0*field(i,j-0,k)
          + field(i,j+1,k)
        )/dx[1]/dx[1];
      break;
    case 3:
      return (
          field(i,j,k-1)
          - 2.0*field(i,j,k-0)
          + field(i,j,k+1)
        )/dx[2]/dx[2];
      break;
  }

  /* XXX */
  return 0;
}

inline real_t forward_double_derivative_stencil_Odx2(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
          2.0*field(i,j,k)
          - 5.0*field(i+1,j,k)
          + 4.0*field(i+2,j,k)
          - 1.0*field(i+3,j,k)
        )/dx[0]/dx[0];
      break;
    case 2:
      return (
          2.0*field(i,j,k)
          - 5.0*field(i,j+1,k)
          + 4.0*field(i,j+2,k)
          - 1.0*field(i,j+3,k)
        )/dx[1]/dx[1];
      break;
    case 3:
      return (
          2.0*field(i,j,k)
          - 5.0*field(i,j,k+1)
          + 4.0*field(i,j,k+2)
          - 1.0*field(i,j,k+3)
        )/dx[2]/dx[2];
      break;
  }

  /* XXX */
  return 0;
}


inline real_t backward_double_derivative_stencil_Odx2(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
          2.0*field(i,j,k)
          - 5.0*field(i-1,j,k)
          + 4.0*field(i-2,j,k)
          - 1.0*field(i-3,j,k)
        )/dx[0]/dx[0];
      break;
    case 2:
      return (
          2.0*field(i,j,k)
          - 5.0*field(i,j-1,k)
          + 4.0*field(i,j-2,k)
          - 1.0*field(i,j-3,k)
        )/dx[1]/dx[1];
      break;
    case 3:
      return (
          2.0*field(i,j,k)
          - 5.0*field(i,j,k-1)
          + 4.0*field(i,j,k-2)
          - 1.0*field(i,j,k-3)
        )/dx[2]/dx[2];
      break;
  }

  /* XXX */
  return 0;
}

 
inline real_t double_derivative_stencil_Odx4(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
          - 1.0/12.0*field(i-2,j,k)
          + 4.0/3.0*field(i-1,j,k)
          - 5.0/2.0*field(i-0,j,k)
          + 4.0/3.0*field(i+1,j,k)
          - 1.0/12.0*field(i+2,j,k)
        )/dx[0]/dx[0];
      break;
    case 2:
      return (
          - 1.0/12.0*field(i,j-2,k)
          + 4.0/3.0*field(i,j-1,k)
          - 5.0/2.0*field(i,j-0,k)
          + 4.0/3.0*field(i,j+1,k)
          - 1.0/12.0*field(i,j+2,k)
        )/dx[1]/dx[1];
      break;
    case 3:
      return (
          - 1.0/12.0*field(i,j,k-2)
          + 4.0/3.0*field(i,j,k-1)
          - 5.0/2.0*field(i,j,k-0)
          + 4.0/3.0*field(i,j,k+1)
          - 1.0/12.0*field(i,j,k+2)
        )/dx[2]/dx[2];
      break;
  }

  /* XXX */
  return 0;
}

inline real_t double_derivative_stencil_Odx6(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
          1.0/90.0*field(i-3,j,k)
          - 3.0/20.0*field(i-2,j,k)
          + 3.0/2.0*field(i-1,j,k)
          - 49.0/18.0*field(i-0,j,k)
          + 3.0/2.0*field(i+1,j,k)
          - 3.0/20.0*field(i+2,j,k)
          + 1.0/90.0*field(i+3,j,k)
        )/dx[0]/dx[0];
      break;
    case 2:
      return (
          1.0/90.0*field(i,j-3,k)
          - 3.0/20.0*field(i,j-2,k)
          + 3.0/2.0*field(i,j-1,k)
          - 49.0/18.0*field(i,j-0,k)
          + 3.0/2.0*field(i,j+1,k)
          - 3.0/20.0*field(i,j+2,k)
          + 1.0/90.0*field(i,j+3,k)
        )/dx[1]/dx[1];
      break;
    case 3:
      return (
          1.0/90.0*field(i,j,k-3)
          - 3.0/20.0*field(i,j,k-2)
          + 3.0/2.0*field(i,j,k-1)
          - 49.0/18.0*field(i,j,k-0)
          + 3.0/2.0*field(i,j,k+1)
          - 3.0/20.0*field(i,j,k+2)
          + 1.0/90.0*field(i,j,k+3)
        )/dx[2]/dx[2];
      break;
  }

  /* XXX */
  return 0;
}

inline real_t double_derivative_stencil_Odx8(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return (
          - 1.0/560.0*field(i-4,j,k)
          + 8.0/315.0*field(i-3,j,k)
          - 1.0/5.0*field(i-2,j,k)
          + 8.0/5.0*field(i-1,j,k)
          - 205.0/72.0*field(i-0,j,k)
          + 8.0/5.0*field(i+1,j,k)
          - 1.0/5.0*field(i+2,j,k)
          + 8.0/315.0*field(i+3,j,k)
          - 1.0/560.0*field(i+4,j,k)
        )/dx[0]/dx[0];
      break;
    case 2:
      return (
          - 1.0/560.0*field(i,j-4,k)
          + 8.0/315.0*field(i,j-3,k)
          - 1.0/5.0*field(i,j-2,k)
          + 8.0/5.0*field(i,j-1,k)
          - 205.0/72.0*field(i,j-0,k)
          + 8.0/5.0*field(i,j+1,k)
          - 1.0/5.0*field(i,j+2,k)
          + 8.0/315.0*field(i,j+3,k)
          - 1.0/560.0*field(i,j+4,k)
        )/dx[1]/dx[1];
      break;
    case 3:
      return (
          - 1.0/560.0*field(i,j,k-4)
          + 8.0/315.0*field(i,j,k-3)
          - 1.0/5.0*field(i,j,k-2)
          + 8.0/5.0*field(i,j,k-1)
          - 205.0/72.0*field(i,j,k-0)
          + 8.0/5.0*field(i,j,k+1)
          - 1.0/5.0*field(i,j,k+2)
          + 8.0/315.0*field(i,j,k+3)
          - 1.0/560.0*field(i,j,k+4)
        )/dx[2]/dx[2];
      break;
  }

  /* XXX */
  return 0;
}

inline real_t forward_dissipation_stencil_Odx2(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return -1.0/2.0*dx[0]*
        forward_double_derivative_stencil_Odx2(i,j,k,d,field, dx);
      break;
    case 2:
      return -1.0/2.0*dx[1]*
        forward_double_derivative_stencil_Odx2(i,j,k,d,field, dx);
      break;
    case 3:
      return -1.0/2.0*dx[2]*
        forward_double_derivative_stencil_Odx2(i,j,k,d,field, dx);
      break;
  }

  /* XXX */
  return 0;
}


inline real_t backward_dissipation_stencil_Odx2(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  switch (d) {
    case 1:
      return +1.0/2.0*dx[0]*
        backward_double_derivative_stencil_Odx2(i,j,k,d,field, dx);
      break;
    case 2:
      return +1.0/2.0*dx[1]*
        backward_double_derivative_stencil_Odx2(i,j,k,d,field, dx);
      break;
    case 3:
      return +1.0/2.0*dx[2]*
        backward_double_derivative_stencil_Odx2(i,j,k,d,field, dx);
      break;
  }

  /* XXX */
  return 0;
}

inline real_t forward_dissipation_stencil_Odx4(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  /* XXX */
  return 0;
}


inline real_t backward_dissipation_stencil_Odx4(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  /* XXX */
  return 0;
}


inline real_t forward_dissipation_stencil_Odx6(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  /* XXX */
  return 0;
}


inline real_t backward_dissipation_stencil_Odx6(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  /* XXX */
  return 0;
}


inline real_t forward_dissipation_stencil_Odx8(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  /* XXX */
  return 0;
}


inline real_t backward_dissipation_stencil_Odx8(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  /* XXX */
  return 0;
}

 
inline real_t forward_dissipation(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  return STENCIL_ORDER_FUNCTION(forward_dissipation_stencil_Odx)(i, j, k, d, field, dx);
}


inline real_t backward_dissipation(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  return STENCIL_ORDER_FUNCTION(backward_dissipation_stencil_Odx)(i, j, k, d, field, dx);
}



/**
 * @brief Compute a derivative using a stencil order defined by a
 * preprocessor directive
 * 
 * @param i x-index
 * @param j x-index
 * @param k x-index
 * @param d direction of derivative
 * @param field field to differentiate
 * @return derivative
 */
inline real_t derivative(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  return STENCIL_ORDER_FUNCTION(derivative_Odx)(i, j, k, d, field, dx);
}

inline real_t forward_derivative(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  return STENCIL_ORDER_FUNCTION(forward_derivative_Odx)(i, j, k, d, field, dx)
    + STENCIL_ORDER_FUNCTION(forward_dissipation_stencil_Odx)(i, j, k, d, field, dx);
}

inline real_t backward_derivative(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  return STENCIL_ORDER_FUNCTION(backward_derivative_Odx)(i, j, k, d, field, dx)
    + STENCIL_ORDER_FUNCTION(backward_dissipation_stencil_Odx)(i, j, k, d, field, dx);
}

inline real_t upwind_derivative(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[], real_t c)
{
  if( c > 0) return c * forward_derivative(i,j,k,d,field, dx);
  else return c * backward_derivative(i,j,k,d,field, dx);
}

inline real_t bd_derivative(
  idx_t i, idx_t j, idx_t k, int d,
  arr_t & field, const double dx[], idx_t l_idx, idx_t codim)
{
  if(codim == 1)
  {
    switch (d) {
    case 1:
      if(l_idx==0  )
        return forward_derivative(i,j,k,d,field,dx);
      else if(l_idx == 1)
        return backward_derivative(i,j,k,d,field,dx);
      else
        return derivative(i,j,k,d,field,dx);
      break;
    case 2:
      if(l_idx == 2)
        return forward_derivative(i,j,k,d,field, dx);
      else if(l_idx == 3)
        return backward_derivative(i,j,k,d,field, dx);
      else
        return derivative(i,j,k,d,field,dx);
      break;
    case 3:
      if( l_idx == 4)
        return forward_derivative(i,j,k,d,field, dx);
      else if(l_idx == 5)
        return backward_derivative(i,j,k,d,field, dx);
      else
        return derivative(i,j,k,d,field,dx);
      break;
    }
  }
  else if(codim == 2)
  {
    switch (d) {
    case 1:
      if(l_idx == 0 || l_idx == 2 || l_idx == 4 || l_idx == 6)
        return forward_derivative(i,j,k,d,field,dx);
      else if(l_idx == 1 || l_idx == 3 || l_idx == 5 || l_idx == 7)
        return backward_derivative(i,j,k,d,field,dx);
      else
        return derivative(i,j,k,d,field,dx);
      break;
    case 2:
      if(l_idx == 0 || l_idx == 1 || l_idx == 8 || l_idx == 10)
        return forward_derivative(i,j,k,d,field, dx);
      else if(l_idx == 2 || l_idx == 3 || l_idx ==9 || l_idx == 11)
        return backward_derivative(i,j,k,d,field, dx);
      else
        return derivative(i,j,k,d,field,dx);
      break;
    case 3:
      if( l_idx == 4 || l_idx == 5 || l_idx == 8 || l_idx == 9)
        return forward_derivative(i,j,k,d,field, dx);
      else if(l_idx == 6 || l_idx == 7 || l_idx == 10 || l_idx == 11)
        return backward_derivative(i,j,k,d,field, dx);
      else
        return derivative(i,j,k,d,field,dx);
      break;
    }
    
  }
  else if(codim == 3)
  {
    switch (d) {
    case 1:
      if(l_idx == 0 || l_idx == 2 || l_idx == 4 || l_idx == 6)
        return forward_derivative(i,j,k,d,field,dx);
      else if(l_idx == 1 || l_idx == 3 || l_idx == 5 || l_idx == 7)
        return backward_derivative(i,j,k,d,field,dx);
      else
        return derivative(i,j,k,d,field,dx);
      break;
    case 2:
      if(l_idx == 0 || l_idx == 1 || l_idx == 4 || l_idx == 5)
        return forward_derivative(i,j,k,d,field, dx);
      else if(l_idx == 2 || l_idx == 3 || l_idx ==6 || l_idx == 7)
        return backward_derivative(i,j,k,d,field, dx);
      else
        return derivative(i,j,k,d,field,dx);
      break;
    case 3:
      if( l_idx == 0 || l_idx == 1 || l_idx == 2 || l_idx == 3)
        return forward_derivative(i,j,k,d,field, dx);
      else if(l_idx == 4 || l_idx == 5 || l_idx == 6 || l_idx == 7)
        return backward_derivative(i,j,k,d,field, dx);
      else
        return derivative(i,j,k,d,field,dx);
      break;
    }
    
  }
  else
    TBOX_ERROR("Codim value is not in the range of 1-3\n");

  /* XXX */
  return 0;
}

inline real_t derivative_norm(idx_t i, idx_t j, idx_t k,
  arr_t & field)
{
  return sqrt(pw2(field(i+1,j,k) - field(i-1,j,k))
             + pw2(field(i,j+1,k) - field(i,j-1,k))
              + pw2(field(i,j,k+1) - field(i,j,k-1)));
}
 
/**
 * @brief Compute a mixed derivative using a stencil order defined by a
 * preprocessor directive
 * 
 * @param i x-index
 * @param j x-index
 * @param k x-index
 * @param d1 direction of derivative in one direction
 * @param d2 direction of derivative in another direction
 * @param field field to differentiate
 * @return derivative
 */
inline real_t mixed_derivative_stencil(idx_t i, idx_t j, idx_t k, int d1, int d2, arr_t & field, const double dx[])
{
  return STENCIL_ORDER_FUNCTION(mixed_derivative_stencil_Odx)(i, j, k, d1, d2, field, dx);
}

/**
 * @brief Compute a second-order derivative using a stencil order defined by a
 * preprocessor directive
 * 
 * @param i x-index
 * @param j x-index
 * @param k x-index
 * @param d direction of 2nd order derivative to compute
 * @param field field to differentiate
 * @return derivative
 */
inline real_t double_derivative_stencil(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field, const double dx[])
{
  return STENCIL_ORDER_FUNCTION(double_derivative_stencil_Odx)(i, j, k, d, field, dx);
}

/**
 * @brief A more generic function for 2nd derivs; calls
 * either @double_derivative_stencil or @mixed_derivative_stencil
 * 
 * @param i x-index
 * @param j x-index
 * @param k x-index
 * @param d1 direction of first derivative
 * @param d2 direction of second derivative
 * @param field field to differentiate
 * @return derivative
 */
inline real_t double_derivative(idx_t i, idx_t j, idx_t k, int d1, int d2,
    arr_t & field, const double dx[])
{
  if(d1 == d2) {
    return double_derivative_stencil(i, j, k, d1, field, dx);
  } else {
    return mixed_derivative_stencil(i, j, k, d1, d2, field, dx);
  }

  /* XXX */
  return 0;
}

/**
 * @brief Computes the laplacian of a field
 * @details sums up double_derivatives
 * 
 * @brief A more generic function for 2nd derivs; calls
 * either @double_derivative_stencil or @mixed_derivative_stencil
 * 
 * @param i x-index
 * @param j x-index
 * @param k x-index
 * @return laplacian
 */
inline real_t laplacian(idx_t i, idx_t j, idx_t k, arr_t & field, const double dx[])
{
  return (
    double_derivative(i, j, k, 1, 1, field, dx)
    + double_derivative(i, j, k, 2, 2, field, dx)
    + double_derivative(i, j, k, 3, 3, field, dx)
  );
}

/**
 * @brief Compute the average of a field
 * 
 * @param field field to average
 * @return average
 */
/* inline real_t average(arr_t & field) */
/* { */
/*   // note this may have poor precision for large datasets */
/*   real_t sum = 0.0;  */
/*   idx_t i=0, j=0, k=0; */
  
/*   LOOP3(i, j, k) */
/*   { */
/*     sum += field(i,j,k); */
/*   } */
/*   return sum/POINTS; */
/* } */

/* /\** */
/*  * @brief Compute the average volume of a simulation */
/*  *  */
/*  * @param DIFFphi difference variable conformal factor field */
/*  * @param phi_FRW reference FRW conformal factor */
/*  *  */
/*  * @return [description] */
/*  *\/ */
/* inline real_t volume_average(arr_t & DIFFphi, real_t phi_FRW) */
/* { */
/*   real_t sum = 0.0;  */
/*   idx_t i=0, j=0, k=0; */

/*   #pragma omp parallel for default(shared) private(i, j, k) reduction(+:sum) */
/*   LOOP3(i, j, k) */
/*   { */
/*     sum += exp(6.0*(DIFFphi(i,j,k) + phi_FRW)); */
/*   } */
/*   return sum/POINTS; */
/* } */

/* /\** */
/*  * @brief Compute the volume-weighted average of a field */
/*  *  */
/*  * @param field field to average */
/*  * @param DIFFphi difference variable conformal factor field */
/*  * @param phi_FRW reference FRW conformal factor */
/*  * @return volume-weighted average */
/*  *\/ */
/* inline real_t conformal_average(arr_t & field, arr_t & DIFFphi, real_t phi_FRW) */
/* { */
/*   real_t sum = 0.0;  */
/*   idx_t i=0, j=0, k=0; */

/*   #pragma omp parallel for default(shared) private(i, j, k) reduction(+:sum) */
/*   LOOP3(i, j, k) */
/*   { */
/*     sum += exp(6.0*(DIFFphi(i,j,k) + phi_FRW))*field(i,j,k); */
/*   } */
/*   real_t vol = volume_average(DIFFphi, phi_FRW); */
/*   return sum/POINTS/vol; */
/* } */

/* /\** */
/*  * @brief Compute the standard deviation of a field */
/*  *  */
/*  * @param field field to analyze */
/*  * @param avg pre-computed average of a field */
/*  * @return standard deviation */
/*  *\/ */
/* inline real_t standard_deviation(arr_t & field, real_t avg) */
/* { */
/*   // note this may have poor precision for large datasets */
/*   idx_t i=0, j=0, k=0; */
/*   real_t sum = 0.0; */
/*   #pragma omp parallel for default(shared) private(i, j, k) reduction(+:sum) */
/*   LOOP3(i, j, k) */
/*   { */
/*     sum += pw2(avg - field(i,j,k)); */
/*   } */
/*   return sqrt(sum/(POINTS-1)); */
/* } */

/* /\** */
/*  * @brief Compute the standard deviation of a field */
/*  *  */
/*  * @param field field to analyze */
/*  * @return standard deviation */
/*  *\/ */
/* inline real_t standard_deviation(arr_t & field) */
/* { */
/*   real_t avg = average(field); */
/*   return standard_deviation(field, avg); */
/* } */

/* /\** */
/*  * @brief Compute the volume-weighted standard deviation of a field */
/*  *  */
/*  * @param field field to analyze */
/*  * @param DIFFphi difference variable conformal factor field */
/*  * @param phi_FRW reference FRW conformal factor */
/*  * @param avg pre-computed volume-weighted average */
/*  * @return volume-weighted standard deviation */
/*  *\/ */
/* inline real_t conformal_standard_deviation( */
/*   arr_t & field, arr_t & DIFFphi, real_t phi_FRW, real_t avg) */
/* { */
/*   real_t sum = 0.0;  */
/*   idx_t i=0, j=0, k=0; */

/*   #pragma omp parallel for default(shared) private(i, j, k) reduction(+:sum) */
/*   LOOP3(i, j, k) */
/*   { */
/*     sum += exp(6.0*(DIFFphi(i,j,k) + phi_FRW))*pw2(avg - field(i,j,k)); */
/*   } */
/*   real_t vol = volume_average(DIFFphi, phi_FRW); */
/*   return sqrt(sum/(POINTS-1)/vol); */
/* } */

/* /\** */
/*  * @brief Compute the volume-weighted standard deviation of a field */
/*  *  */
/*  * @param field field to analyze */
/*  * @param DIFFphi difference variable conformal factor field */
/*  * @param phi_FRW reference FRW conformal factor */
/*  * @return volume-weighted standard deviation */
/*  *\/ */
/* inline real_t conformal_standard_deviation( */
/*   arr_t & field, arr_t & DIFFphi, real_t phi_FRW) */
/* { */
/*   real_t avg = conformal_average(field, DIFFphi, phi_FRW); */
/*   return conformal_standard_deviation(field, DIFFphi, phi_FRW, avg); */
/* } */

/* /\** */
/*  * @brief Compute the maximum value of a field */
/*  *\/ */
/* inline real_t max(arr_t & field) */
/* { */
/*   idx_t i=0, j=0, k=0; */
/*   real_t max_val = field[0]; */
/*   LOOP3(i, j, k) */
/*   { */
/*     if(field(i,j,k) > max_val) { */
/*       max_val = field(i,j,k); */
/*     } */
/*   } */
/*   return max_val; */
/* } */

/* /\** */
/*  * @brief Compute the minimum value of a field */
/*  *\/ */
/* inline real_t min(arr_t & field) */
/* { */
/*   idx_t i=0, j=0, k=0; */
/*   real_t min_val = field[0]; */
/*   LOOP3(i, j, k) */
/*   { */
/*     if(field(i,j,k) < min_val) { */
/*       min_val = field(i,j,k); */
/*     } */
/*   } */
/*   return min_val; */
/* } */

/* /\** */
/*  * @brief compute the number of NAN values in a field */
/*  * @details may have portability issues */
/*  *\/ */
/* inline idx_t numNaNs(arr_t & field) */
/* { */
/*   idx_t i=0, j=0, k=0; */
/*   idx_t NaNs = 0; */
  
/*   LOOP3(i,j,k) */
/*   { */
/*     real_t val = field(i,j,k); */
/*     union { float val; uint32_t x; } u = { (float) val }; */
/*     if((u.x << 1) > 0xff000000u) */
/*     { */
/*       NaNs += 1; */
/*     } */
/*   } */

/*   return NaNs; */
/* } */


/* inline idx_t idx_t_mod(idx_t n, idx_t d) */
/* { */
/*   idx_t mod = n % d; */
/*   if(mod < 0) */
/*     mod += d; */
/*   return mod; */
/* } */

/* inline real_t real_t_mod(real_t n, real_t d) */
/* { */
/*   return n - d*std::floor(n/d); */
/* } */





}
#endif
