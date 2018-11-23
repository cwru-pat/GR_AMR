#ifndef FAS_MULTIGRID_H
#define FAS_MULTIGRID_H

#include <omp.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include "../../utils/Array.h"

#include "multigrid_bd_handler.h"

#include "../../cosmo_macros.h"

#define PI  (4.0*atan(1.0))

#define FAS_LOOP3_N(i, j, k, nx, ny, nz)  \
  for(i=0; i<nx; ++i)                     \
    for(j=0; j<ny; ++j)                   \
      for(k=0; k<nz; ++k)



using namespace SAMRAI;

namespace cosmo
{

/**
 * @brief single element in a term
 */
typedef struct{
  idx_t type;    ///< element type; 0 for constant function, 1 for polynomial, 2-10 for single and double derivatives, 11 for laplacian (see also enum atom_type)
  idx_t u_id;    ///< id of varible needs to be solved, won't be visited when type = 0
  real_t value;  ///< exponent value, has meaning only when type = 1 (polynomial)
} atom;


/**
 * @brief single term in a differential equation
 * @details
 * combined by multiplication of "atoms"
 */
class molecule
{
 public:
  atom * atoms;      ///< vector storing "atoms"
  idx_t atom_n;      ///< "atom" number
  real_t const_coef; ///< constant coefficient of single term

  molecule()
  {
    atom_n = 0;
  }

  ~molecule()
  {
    delete [] atoms;
  }

  void init(idx_t atom_n_in, real_t const_coef_in)
  {
    atom_n = 0;
    atoms = new atom[atom_n_in];

    const_coef = const_coef_in;
  }

  void add_atom(atom atom_in)
  {
    atoms[atom_n++] = atom_in;
  }
};

class FASMultigrid
{
  private:

  // grid (array) type
  typedef CosmoArray<idx_t, real_t> fas_grid_t;
  // heirarchy type (set of some grids at different depths)
  typedef fas_grid_t * fas_heirarchy_t;
  // set of heirarchies (one for each variable/equation)
  typedef fas_heirarchy_t * fas_heirarchy_set_t;

  // define heirarchy of references to grids
  fas_heirarchy_set_t u_h;             ///< field seeking a solution for
  fas_heirarchy_set_t tmp_h;           ///< reusable grid for storing intermediate calculations
  fas_heirarchy_set_t coarse_src_h;    ///< multigrid source term
  fas_heirarchy_set_t jac_rhs_h;       ///< - F(u) which is rhs of Jacob Linear function
  fas_heirarchy_set_t damping_v_h;     ///< _lap (u) - f, used to calculate F(u + \lambda v)
  fas_heirarchy_set_t * rho_h;         ///< source matter terms with number being rho_num;
  fas_heirarchy_set_t u_shift_h;    ///< the shift function of u, only change vars in poly to (u+u_shift) 
  
  idx_t u_n;          ///< number of variables ( = number of equations)

  idx_t * molecule_n; ///< number of molecules for each equation

  idx_t *nx_h, *ny_h, *nz_h;  ///< number of grid points in each direction at different depths

  real_t relaxation_tolerance;  ///< desired precision when performing relaxation

  idx_t max_depth, max_depth_idx;
  idx_t min_depth, min_depth_idx;
  idx_t total_depths, max_relax_iters;

  idx_t der_type[12 /* number of items in enum atom_type */][2 /* derivative directions(s) */];      ///< vectors that stores devivative directions

  real_t double_der_coef[9];  ///< vectors that stores coefficients of f(x,y,z) for different order stencils, used for jac equation iteration

  real_t H_LEN_FRAC[3];

  std::string boundary_type;

  multigridBdHandler * bd_handler;
    
  /**
   * @brief indexing scheme of a grid heirarchy
   * @description return index of grid at a particular depth
   *  in a grid heirarchy
   * 
   * @param depth "depth" of grid
   * @return index
   */
  inline idx_t _dIdx(idx_t depth)
  {
    return depth - min_depth;
  }

  /**
   * @brief return sign of argument
   * @details return zerp when argument is zero
   */
  inline idx_t _sign(real_t x)
  {
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
  }
  
  /**
   * @brief compute power of 2
   * 
   * @param pwr power to raise 2 to
   * @return 2^pwr
   */
  inline idx_t _2toPwr(idx_t pwr)
  {
    return 1<<pwr;
  }

  /**
   * @brief compute integer number to power of 3
   * 
   * @param number to raise to ^3
   * @return pwr^3
   */
  inline idx_t _Pwr3(idx_t num)
  {
    return num*num*num;
  }

  /**
   * @brief compute integer number to power of 2
   * 
   * @param number to raise to ^2
   * @return pwr^2
   */
  inline real_t _Pwr2(real_t num)
  {
    return num * num;
  }

 public:
  
  // enum for relaxation type
  enum relax_t
  {
    inexact_newton,
    inexact_newton_constrained, // inexact Newton with volume constraint enforced
    newton
  };

  relax_t relax_scheme;
  
  enum atom_type
  {
    poly = 1,
    der1 = 2,
    der2 = 3,
    der3 = 4,
    der11 = 5,
    der22 = 6,
    der33 = 7,
    der12 = 8,
    der13 = 9,
    der23 = 10,
    lap = 11
  };

  molecule ** eqns; ///< All terms in all equations

  char * interpolator_pars;

  int interp_param_handle;
  
  FASMultigrid(fas_grid_t u_in[], idx_t u_n_in, idx_t molecule_n_in [],
               idx_t max_depth_in, idx_t max_relax_iters_in,
               real_t relaxation_tolerance_in, real_t H_LEN_FRAC_IN[3], idx_t NX, idx_t NY, idx_t NZ, multigridBdHandler * bd_handler_in);
  ~FASMultigrid();

  void add_atom_to_eqn(atom atom_in, idx_t molecule_id, idx_t eqn_id);

  real_t _evaluateEllipticEquationPt(idx_t eqn_id, idx_t depth_idx, idx_t i,
    idx_t j, idx_t k);

  void _evaluateIterationForJacEquation(idx_t eqn_id, idx_t depth_idx,
    real_t &coef_a, real_t &coef_b, idx_t i, idx_t j, idx_t k, idx_t u_id);

  real_t _evaluateDerEllipticEquation(idx_t eqn_id, idx_t depth_idx, idx_t i,
    idx_t j, idx_t k, idx_t var_id);

  void _zeroGrid(fas_grid_t & grid);

  real_t _totalGrid(fas_grid_t & grid);

  real_t _averageGrid(fas_grid_t & grid);

  real_t _maxGrid(fas_grid_t &grid);

  real_t _minGrid(fas_grid_t & grid);

  void _shiftGridVals(fas_grid_t & grid, real_t shift);

  void _restrictFine2coarse(fas_heirarchy_t grid_heirarchy, idx_t fine_depth);

  void _interpolateCoarse2fine(fas_heirarchy_t grid_heirarchy,
    idx_t coarse_depth);

  void _evaluateEllipticEquation(fas_heirarchy_t  result_h, idx_t eqn_id,
    idx_t depth);

  void _computeResidual(fas_heirarchy_t residual_h, idx_t eqn_id, idx_t depth);

  real_t _getMaxResidual(idx_t eqn_id, idx_t depth);

  real_t _getMaxResidualAllEqs(idx_t depth);

  void _computeCoarseRestrictions(idx_t eqn_id, idx_t fine_depth);

  void _changeApproximateSolutionToError(fas_heirarchy_t  appx_to_err_h,
    fas_heirarchy_t  exact_soln_h, idx_t depth);

  void _correctFineFromCoarseErr_Err2Appx(fas_heirarchy_t err2appx_h,
    fas_heirarchy_t  appx_soln_h, idx_t fine_depth);

  void _copyGrid(fas_heirarchy_t from_h[], fas_heirarchy_t to_h[],
    idx_t eqn_id, idx_t depth);

  bool _getLambda( idx_t depth, real_t norm);

  bool _jacobianRelax( idx_t depth, real_t norm, real_t C, idx_t p);

  bool _singularityExists(idx_t eqn_id, idx_t depth);

  bool _relaxSolution_GaussSeidel( idx_t depth, idx_t max_iterations);

  void _printStrip(fas_grid_t & out_h);

  void build_rho();

  bool VCycle();

  void VCycles(idx_t num_cycles);

  void setPolySrcAtPt(idx_t eqn_id, idx_t mol_id, idx_t i, idx_t j, idx_t k,
    real_t value);

  void setShiftSrcAtPt(idx_t u_id, idx_t i, idx_t j, idx_t k, real_t value);
  
  void initializeRhoHeirarchy();

  void initializeShiftHeirarchy();
  
  void printSolutionStrip(idx_t depth);

  real_t double_derivative(
    idx_t i, idx_t j, idx_t k, idx_t nx, idx_t ny, idx_t nz,
    int d1, int d2, fas_grid_t & field);

  
  real_t double_derivative_stencil(
    idx_t i, idx_t j, idx_t k, idx_t nx, idx_t ny, idx_t nz, int d, fas_grid_t & field);

  real_t mixed_derivative_stencil(
    idx_t i, idx_t j, idx_t k,idx_t nx, idx_t ny, idx_t nz,
    int d1, int d2, fas_grid_t & field);

  
  real_t laplacian(
    idx_t i, idx_t j, idx_t k, idx_t nx, idx_t ny, idx_t nz, fas_grid_t & field);

  real_t derivative(
    idx_t i, idx_t j, idx_t k, idx_t nx, idx_t ny, idx_t nz, int d, fas_grid_t & field);

  void fillBoundary(fas_grid_t & data);
};

} // namespace cosmo
#endif
