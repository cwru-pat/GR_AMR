#ifndef COSMO_ICS
#define COSMO_ICS

#define NP_INDEX(i,j,k) ((NZ)*(NY)*(i) + (NZ)*(j) + (k))
#define FFT_NP_INDEX(i,j,k) ((NZ/2+1)*NY*(i) + (NZ/2+1)*(j) + (k))

#include "../cosmo_includes.h"
#include "ICs_data.h"
#include <fftw3.h>
#include <zlib.h>

using namespace SAMRAI;

namespace cosmo
{

ICsData cosmo_get_ICsData(
  std::shared_ptr<tbox::Database> cosmo_ICs_db, real_t domain_size);

real_t cosmo_power_spectrum(real_t k, ICsData *icd);
void set_gaussian_random_field(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, idx_t f_id, ICsData *icd);
 

} // namespace cosmo

#endif
