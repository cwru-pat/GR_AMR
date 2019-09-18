#ifndef COSMO_DUST_FLUID_DATA
#define COSMO_DUST_FLUID_DATA

namespace cosmo
{
typedef struct {

  real_t D, S1, S2, S3, E;
  real_t Si1, Si2, Si3;
  real_t v1, v2, v3, vi1, vi2, vi3;
  real_t m11, m22, m33, m12, m13, m23;
  real_t d1m11, d1m22, d1m33, d1m12, d1m13, d1m23;
  real_t d2m11, d2m22, d2m33, d2m12, d2m13, d2m23;
  real_t d3m11, d3m22, d3m33, d3m12, d3m13, d3m23;
  real_t S11, S22, S33, S12, S13, S23;
  real_t K11, K22, K33, K12, K13, K23;
} DustFluidData;

}
#endif
