#ifndef COSMO_GEODESIC_DATA
#define COSMO_GEODESIC_DATA

namespace cosmo
{

typedef struct {
  real_t m11, m12, m13, m22, m23, m33;
  real_t mi11, mi12, mi13, mi22, mi23, mi33;
  real_t q1, q2, q3, qi1, qi2, qi3;

  real_t chi;
  real_t d1chi, d2chi, d3chi;

  real_t beta1, beta2, beta3;
  real_t d1beta1, d1beta2, d1beta3;
  real_t d2beta1, d2beta2, d2beta3;
  real_t d3beta1, d3beta2, d3beta3;

  real_t alpha, d1alpha, d2alpha, d3alpha;

  real_t d1m11, d1m12, d1m13, d1m22, d1m23, d1m33;
  real_t d2m11, d2m12, d2m13, d2m22, d2m23, d2m33;
  real_t d3m11, d3m12, d3m13, d3m22, d3m23, d3m33;

  real_t p0;
  double x, y, z;
  double lambda;
  
}GeodesicData;

 
}
#endif
