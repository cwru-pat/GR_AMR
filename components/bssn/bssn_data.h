#ifndef COSMO_BSSN_DATA
#define COSMO_BSSN_DATA

#include "../../cosmo_includes.h"
#include "../../cosmo_types.h"
#include "bssn_macros.h"

namespace cosmo
{

/**
 * @struct BSSNData
 * @brief Structure containing BSSN metric variables and various derived
 * quantities, such as derivatives of BSSN variables, christoffel symbols,
 * etc. Most undocumented variables correspond to values taken from a
 * particular field; most derived variables are documented.
 */
typedef struct {

  idx_t i, j, k, idx;

  // local copies of current field values
  BSSN_APPLY_TO_FIELDS(DECLARE_REAL_T)
  // Source terms
  BSSN_APPLY_TO_SOURCES(DECLARE_REAL_T)
  // "extra" fields
  BSSN_APPLY_TO_GEN1_EXTRAS(DECLARE_REAL_T)

  // non-differenced quantities
  real_t chi; ///< conformal factor \f$\phi\f$
  real_t K; ///< extrinsic curvature \f$K\f$
  real_t r; ///< density \f$\rho\f$
  real_t S; ///< "pressure" (trace of \f$S_{ij}\f$, or \f$\gamma^{ij} S_{ij}\f$)
  real_t alpha; ///< lapse, \f$\alpha\f$
  real_t gamma11, ///< \f$\bar{\gamma}_{11}\f$ (conformal 11 metric component)
         gamma12, ///< \f$\bar{\gamma}_{12}\f$ (conformal 12 metric component)
         gamma13, ///< \f$\bar{\gamma}_{13}\f$ (conformal 13 metric component)
         gamma22, ///< \f$\bar{\gamma}_{22}\f$ (conformal 22 metric component)
         gamma23, ///< \f$\bar{\gamma}_{23}\f$ (conformal 23 metric component)
         gamma33; ///< \f$\bar{\gamma}_{33}\f$ (conformal 33 metric component)
  real_t gammai11, ///< \f$\bar{\gamma}^{11}\f$ (inverse conformal 11 metric component)
         gammai12, ///< \f$\bar{\gamma}^{12}\f$ (inverse conformal 12 metric component)
         gammai13, ///< \f$\bar{\gamma}^{13}\f$ (inverse conformal 13 metric component)
         gammai22, ///< \f$\bar{\gamma}^{22}\f$ (inverse conformal 22 metric component)
         gammai23, ///< \f$\bar{\gamma}^{23}\f$ (inverse conformal 23 metric component)
         gammai33; ///< \f$\bar{\gamma}^{33}\f$ (inverse conformal 33 metric component)

  // generic var for misc. expressions
  real_t trace, ///< Generic re-usable variable for trace
         expression; ///< Generic re-usable variable for any expression

  // ricci tensor components
  real_t ricci11, ///< full Ricci tensor component, \f$R_{11}\f$
         ricci12, ///< full Ricci tensor component, \f$R_{12}\f$
         ricci13, ///< full Ricci tensor component, \f$R_{13}\f$
         ricci22, ///< full Ricci tensor component, \f$R_{22}\f$
         ricci23, ///< full Ricci tensor component, \f$R_{23}\f$
         ricci33; ///< full Ricci tensor component, \f$R_{33}\f$
  real_t ricciTF11, ///< full trace-free Ricci tensor component, \f$R_{11}^{TF}\f$
         ricciTF12, ///< full trace-free Ricci tensor component, \f$R_{12}^{TF}\f$
         ricciTF13, ///< full trace-free Ricci tensor component, \f$R_{13}^{TF}\f$
         ricciTF22, ///< full trace-free Ricci tensor component, \f$R_{22}^{TF}\f$
         ricciTF23, ///< full trace-free Ricci tensor component, \f$R_{23}^{TF}\f$
         ricciTF33; ///< full trace-free Ricci tensor component, \f$R_{33}^{TF}\f$
  real_t Uricci11, ///< unitary Ricci tensor component, \f$\bar{R}_{11}\f$
         Uricci12, ///< unitary Ricci tensor component, \f$\bar{R}_{12}\f$
         Uricci13, ///< unitary Ricci tensor component, \f$\bar{R}_{13}\f$
         Uricci22, ///< unitary Ricci tensor component, \f$\bar{R}_{22}\f$
         Uricci23, ///< unitary Ricci tensor component, \f$\bar{R}_{23}\f$
         Uricci33; ///< unitary Ricci tensor component, \f$\bar{R}_{33}\f$
  real_t unitRicci; ///< unitary Ricci scalar, \f$\bar{R}\f$
  // for constraint checking

  // derivatives of \alpha
    // covariant double-derivatives 
    real_t D1D1aTF, ///< trace-free second covariant derivative of lapse, \f$(D_1 D_1 \alpha)^{TF}\f$
           D1D2aTF, ///< trace-free second covariant derivative of lapse, \f$(D_1 D_2 \alpha)^{TF}\f$
           D1D3aTF, ///< trace-free second covariant derivative of lapse, \f$(D_1 D_3 \alpha)^{TF}\f$
           D2D2aTF, ///< trace-free second covariant derivative of lapse, \f$(D_2 D_2 \alpha)^{TF}\f$
           D2D3aTF, ///< trace-free second covariant derivative of lapse, \f$(D_2 D_3 \alpha)^{TF}\f$
           D3D3aTF; ///< trace-free second covariant derivative of lapse, \f$(D_3 D_3 \alpha)^{TF}\f$
    real_t DDaTR; ///< \f$\gamma^{ij}D_i D_j \alpha\f$
    // normal derivatives of
    real_t d1a, ///< partial of alpha, \f$\partial_1 \alpha\f$
           d2a, ///< partial of alpha, \f$\partial_2 \alpha\f$
           d3a; ///< partial of alpha, \f$\partial_3 \alpha\f$


    // covariant double-derivatives 
    real_t D1D1chi, ///< conformal covariant second derivative of chi, \f$\bar{D}_1 \bar{D}_1 \chi\f$
           D1D2chi, ///< conformal covariant second derivative of chi, \f$\bar{D}_1 \bar{D}_2 \chi\f$
           D1D3chi, ///< conformal covariant second derivative of chi, \f$\bar{D}_1 \bar{D}_3 \chi\f$
           D2D2chi, ///< conformal covariant second derivative of chi, \f$\bar{D}_2 \bar{D}_2 \chi\f$
           D2D3chi, ///< conformal covariant second derivative of chi, \f$\bar{D}_2 \bar{D}_3 \chi\f$
           D3D3chi; ///< conformal covariant second derivative of chi, \f$\bar{D}_3 \bar{D}_3 \chi\f$
    // normal derivatives of
  real_t d1chi, ///< \f$\partial_1 \phi \f$
    d2chi, ///< \f$\partial_2 \phi \f$
    d3chi; ///< \f$\partial_3 \phi \f$

  real_t d1d1chi, ///< partial second derivative of phi, \f$\partial_1 \partial_1 \phi\f$
    d1d2chi, ///< partial second derivative of phi, \f$\partial_1 \partial_2 \phi\f$
    d1d3chi, ///< partial second derivative of phi, \f$\partial_1 \partial_3 \phi\f$
    d2d2chi, ///< partial second derivative of phi, \f$\partial_2 \partial_2 \phi\f$
    d2d3chi, ///< partial second derivative of phi, \f$\partial_2 \partial_3 \phi\f$
    d3d3chi; ///< partial second derivative of phi, \f$\partial_3 \partial_3 \phi\f$

  
  // ders of K
  real_t d1K, ///< \f$\partial_1 K\f$
         d2K, ///< \f$\partial_2 K\f$
         d3K; ///< \f$\partial_3 K\f$

  // Contravariant (upstairs index) ext. curvature
  real_t Acont11, ///< Contravariant form of conformal trace-free extrinsic curvature, \f$ \bar{A}^{11} \f$
         Acont12, ///< Contravariant form of conformal trace-free extrinsic curvature, \f$ \bar{A}^{12} \f$
         Acont13, ///< Contravariant form of conformal trace-free extrinsic curvature, \f$ \bar{A}^{13} \f$
         Acont22, ///< Contravariant form of conformal trace-free extrinsic curvature, \f$ \bar{A}^{22} \f$
         Acont23, ///< Contravariant form of conformal trace-free extrinsic curvature, \f$ \bar{A}^{23} \f$
         Acont33; ///< Contravariant form of conformal trace-free extrinsic curvature, \f$ \bar{A}^{33} \f$

  // Christoffel symbols
  real_t G111, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{1}_{11} \f$
         G112, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{1}_{12} \f$
         G113, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{1}_{13} \f$
         G122, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{1}_{22} \f$
         G123, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{1}_{23} \f$
         G133, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{1}_{33} \f$
         G211, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{2}_{11} \f$
         G212, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{2}_{12} \f$
         G213, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{2}_{13} \f$
         G222, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{2}_{22} \f$
         G223, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{2}_{23} \f$
         G233, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{2}_{33} \f$
         G311, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{3}_{11} \f$
         G312, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{3}_{12} \f$
         G313, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{3}_{13} \f$
         G322, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{3}_{22} \f$
         G323, ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{3}_{23} \f$
         G333; ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{3}_{33} \f$

  // Lowered index christoffel symbols
  real_t GL111, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{111} \f$
         GL112, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{112} \f$
         GL113, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{113} \f$
         GL122, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{122} \f$
         GL123, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{123} \f$
         GL133, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{133} \f$
         GL211, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{211} \f$
         GL212, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{212} \f$
         GL213, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{213} \f$
         GL222, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{222} \f$
         GL223, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{223} \f$
         GL233, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{233} \f$
         GL311, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{311} \f$
         GL312, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{312} \f$
         GL313, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{313} \f$
         GL322, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{322} \f$
         GL323, ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{323} \f$
         GL333; ///< Conformal christoffel symbol of the second kind, \f$ \bar{\Gamma}_{333} \f$

  // contraction of christoffel symbols ("Gamma_d" in Z4c)
  real_t Gammad1, ///< Contraction of christoffel symbol (non-dynamical), \f$\bar{\gamma}^{ij} \bar{\Gamma}^{1}_{ij}\f$
         Gammad2, ///< Contraction of christoffel symbol (non-dynamical), \f$\bar{\gamma}^{ij} \bar{\Gamma}^{2}_{ij}\f$
         Gammad3; ///< Contraction of christoffel symbol (non-dynamical), \f$\bar{\gamma}^{ij} \bar{\Gamma}^{3}_{ij}\f$

  // derivatives of the metric, d_i g_jk
  real_t d1g11, ///< First partial derivative of the conformal metric, \f$\partial_1 \bar{\gamma}_{11} \f$
         d1g12, ///< First partial derivative of the conformal metric, \f$\partial_1 \bar{\gamma}_{12} \f$
         d1g13, ///< First partial derivative of the conformal metric, \f$\partial_1 \bar{\gamma}_{13} \f$
         d1g22, ///< First partial derivative of the conformal metric, \f$\partial_1 \bar{\gamma}_{22} \f$
         d1g23, ///< First partial derivative of the conformal metric, \f$\partial_1 \bar{\gamma}_{23} \f$
         d1g33, ///< First partial derivative of the conformal metric, \f$\partial_1 \bar{\gamma}_{33} \f$
         d2g11, ///< First partial derivative of the conformal metric, \f$\partial_2 \bar{\gamma}_{11} \f$
         d2g12, ///< First partial derivative of the conformal metric, \f$\partial_2 \bar{\gamma}_{12} \f$
         d2g13, ///< First partial derivative of the conformal metric, \f$\partial_2 \bar{\gamma}_{13} \f$
         d2g22, ///< First partial derivative of the conformal metric, \f$\partial_2 \bar{\gamma}_{22} \f$
         d2g23, ///< First partial derivative of the conformal metric, \f$\partial_2 \bar{\gamma}_{23} \f$
         d2g33, ///< First partial derivative of the conformal metric, \f$\partial_2 \bar{\gamma}_{33} \f$
         d3g11, ///< First partial derivative of the conformal metric, \f$\partial_3 \bar{\gamma}_{11} \f$
         d3g12, ///< First partial derivative of the conformal metric, \f$\partial_3 \bar{\gamma}_{12} \f$
         d3g13, ///< First partial derivative of the conformal metric, \f$\partial_3 \bar{\gamma}_{13} \f$
         d3g22, ///< First partial derivative of the conformal metric, \f$\partial_3 \bar{\gamma}_{22} \f$
         d3g23, ///< First partial derivative of the conformal metric, \f$\partial_3 \bar{\gamma}_{23} \f$
         d3g33; ///< First partial derivative of the conformal metric, \f$\partial_3 \bar{\gamma}_{33} \f$


  // second derivatives of the metric d_i d_j g_kl
  real_t d1d1g11, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_1 \bar{\gamma}_{11}\f$
         d1d1g12, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_1 \bar{\gamma}_{12}\f$
         d1d1g13, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_1 \bar{\gamma}_{13}\f$
         d1d1g22, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_1 \bar{\gamma}_{22}\f$
         d1d1g23, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_1 \bar{\gamma}_{23}\f$
         d1d1g33, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_1 \bar{\gamma}_{33}\f$
         d1d2g11, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_2 \bar{\gamma}_{11}\f$
         d1d2g12, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_2 \bar{\gamma}_{12}\f$
         d1d2g13, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_2 \bar{\gamma}_{13}\f$
         d1d2g22, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_2 \bar{\gamma}_{22}\f$
         d1d2g23, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_2 \bar{\gamma}_{23}\f$
         d1d2g33, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_2 \bar{\gamma}_{33}\f$
         d1d3g11, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_3 \bar{\gamma}_{11}\f$
         d1d3g12, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_3 \bar{\gamma}_{12}\f$
         d1d3g13, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_3 \bar{\gamma}_{13}\f$
         d1d3g22, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_3 \bar{\gamma}_{22}\f$
         d1d3g23, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_3 \bar{\gamma}_{23}\f$
         d1d3g33, ///< Second partial derivative of the conformal metric, \f$\partial_1 \partial_3 \bar{\gamma}_{33}\f$
         d2d2g11, ///< Second partial derivative of the conformal metric, \f$\partial_2 \partial_2 \bar{\gamma}_{11}\f$
         d2d2g12, ///< Second partial derivative of the conformal metric, \f$\partial_2 \partial_2 \bar{\gamma}_{12}\f$
         d2d2g13, ///< Second partial derivative of the conformal metric, \f$\partial_2 \partial_2 \bar{\gamma}_{13}\f$
         d2d2g22, ///< Second partial derivative of the conformal metric, \f$\partial_2 \partial_2 \bar{\gamma}_{22}\f$
         d2d2g23, ///< Second partial derivative of the conformal metric, \f$\partial_2 \partial_2 \bar{\gamma}_{23}\f$
         d2d2g33, ///< Second partial derivative of the conformal metric, \f$\partial_2 \partial_2 \bar{\gamma}_{33}\f$
         d2d3g11, ///< Second partial derivative of the conformal metric, \f$\partial_2 \partial_3 \bar{\gamma}_{11}\f$
         d2d3g12, ///< Second partial derivative of the conformal metric, \f$\partial_2 \partial_3 \bar{\gamma}_{12}\f$
         d2d3g13, ///< Second partial derivative of the conformal metric, \f$\partial_2 \partial_3 \bar{\gamma}_{13}\f$
         d2d3g22, ///< Second partial derivative of the conformal metric, \f$\partial_2 \partial_3 \bar{\gamma}_{22}\f$
         d2d3g23, ///< Second partial derivative of the conformal metric, \f$\partial_2 \partial_3 \bar{\gamma}_{23}\f$
         d2d3g33, ///< Second partial derivative of the conformal metric, \f$\partial_2 \partial_3 \bar{\gamma}_{33}\f$
         d3d3g11, ///< Second partial derivative of the conformal metric, \f$\partial_3 \partial_3 \bar{\gamma}_{11}\f$
         d3d3g12, ///< Second partial derivative of the conformal metric, \f$\partial_3 \partial_3 \bar{\gamma}_{12}\f$
         d3d3g13, ///< Second partial derivative of the conformal metric, \f$\partial_3 \partial_3 \bar{\gamma}_{13}\f$
         d3d3g22, ///< Second partial derivative of the conformal metric, \f$\partial_3 \partial_3 \bar{\gamma}_{22}\f$
         d3d3g23, ///< Second partial derivative of the conformal metric, \f$\partial_3 \partial_3 \bar{\gamma}_{23}\f$
         d3d3g33; ///< Second partial derivative of the conformal metric, \f$\partial_3 \partial_3 \bar{\gamma}_{33}\f$


  // H constraint calc.
  real_t H; ///< Hamiltonian constraint violation
  // Misc. debugging calc
  real_t db; ///< Misc. re-usable debugging variable

  // additional variables to handle absence of Z4c terms in macros
  // (make sure these get initialized to 0!)
  real_t d1theta, ///< \f$\partial_1 \theta\f$ variable
         d2theta, ///< \f$\partial_2 \theta\f$ variable
         d3theta; ///< \f$\partial_3 \theta\f$ variable

  real_t d1thetau, ///< \f$\partial_1 \theta\f$ variable
         d2thetau, ///< \f$\partial_2 \theta\f$ variable
         d3thetau; ///< \f$\partial_3 \theta\f$ variable

  // additional variables to handle absence of shift terms in macros
  // (make sure these get initialized to 0!)
  real_t d1beta1, ///< derivative of shift, \f$ \partial_1 \beta^1 \f$
         d2beta1, ///< derivative of shift, \f$ \partial_2 \beta^1 \f$
         d3beta1; ///< derivative of shift, \f$ \partial_3 \beta^1 \f$
  real_t d1beta2, ///< derivative of shift, \f$ \partial_1 \beta^2 \f$
         d2beta2, ///< derivative of shift, \f$ \partial_2 \beta^2 \f$
         d3beta2; ///< derivative of shift, \f$ \partial_3 \beta^2 \f$
  real_t d1beta3, ///< derivative of shift, \f$ \partial_1 \beta^3 \f$
         d2beta3, ///< derivative of shift, \f$ \partial_2 \beta^3 \f$
         d3beta3; ///< derivative of shift, \f$ \partial_3 \beta^3 \f$
  #if !USE_BSSN_SHIFT
    real_t beta1; ///< shift, \f$\beta^1\f$
    real_t beta2; ///< shift, \f$\beta^2\f$
    real_t beta3; ///< shift, \f$\beta^3\f$
  #endif

  /* #if USE_EXPANSION */
  /* real_t d1expN, */
  /*        d2expN, */
  /*        d3expN; */
  /* real_t d1expNu, */
  /*        d2expNu, */
  /*        d3expNu; */
  /* #endif */
  // Reference FRW quantities
  real_t chi_FRW; ///< Reference FRW variable, \f$\phi_{FRW}\f$
  real_t K_FRW; ///< Reference FRW variable, \f$K_{FRW}\f$
  real_t rho_FRW, ///< Reference FRW variable, \f$\rho_{FRW}\f$
         S_FRW; ///< Reference FRW variable, \f$S_{FRW}\f$

  real_t K0;
  real_t K_avg;
  
  real_t D1Z1, D2Z2, D3Z3, D1Z2, D1Z3, D2Z3, D2Z1, D3Z1, D3Z2;
  real_t D1Z1TF, D2Z2TF, D3Z3TF, D1Z2TF,
    D1Z3TF, D2Z3TF, D2Z1TF, D3Z1TF, D3Z2TF;
  real_t DZTR;

  real_t norm, x, y, z;

  real_t Z1, Z2, Z3;

  #if USE_DUST_FLUID
  real_t dchidt;
  #endif
  
} BSSNData;

} /* namespace cosmo */

#endif
