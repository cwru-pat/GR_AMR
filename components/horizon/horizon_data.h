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

}
#endif
