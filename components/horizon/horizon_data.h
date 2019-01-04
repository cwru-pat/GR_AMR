#ifndef COSMO_HORIZON_DATA
#define COSMO_HORIZON_DATA

namespace cosmo
{
typedef struct {

  real_t F, s1, s2, s3;
  // field values
  real_t d1F, d2F, d3F;

  real_t d1d1F, d2d2F, d3d3F, d1d2F, d1d3F, d2d3F;
}HorizonData;

typedef struct {

  // Christoffel symbols
  real_t Gc111, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{1}_{11} \f$
         Gc112, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{1}_{12} \f$
         Gc113, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{1}_{13} \f$
         Gc122, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{1}_{22} \f$
         Gc123, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{1}_{23} \f$
         Gc133, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{1}_{33} \f$
         Gc211, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{2}_{11} \f$
         Gc212, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{2}_{12} \f$
         Gc213, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{2}_{13} \f$
         Gc222, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{2}_{22} \f$
         Gc223, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{2}_{23} \f$
         Gc233, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{2}_{33} \f$
         Gc311, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{3}_{11} \f$
         Gc312, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{3}_{12} \f$
         Gc313, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{3}_{13} \f$
         Gc322, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{3}_{22} \f$
         Gc323, ///< Conformal christoffel symbol, \f$ \bar{\Gcamma}^{3}_{23} \f$
         Gc333; ///< Conformal christoffel symbol, \f$ \bar{\Gamma}^{3}_{33} \f$
  
  real_t Gs111,
    Gs112,
    Gs122,
    Gs211,
    Gs212,
    Gs222;

  real_t m11, m12, m13, m22, m23, m33;
  real_t mi11, mi12, mi13, mi22, mi23, mi33;
  real_t q11,
    q12,
    q22;

  real_t qi11,
    qi12,
    qi22;
  
  real_t R;

  real_t R11, R12, R22;

  real_t K11, K12, K13, K22, K23, K33;

  real_t K;

  real_t d1F, d2F, d3F;

  real_t chi;
}KillingData;

 
}
#endif
