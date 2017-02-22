#ifndef COSMO_SCALAR_ICS
#define COSMO_SCALAR_ICS

#include "../../cosmo_includes.h"
#include "scalar.h"
#include <fftw3.h>
#include <zlib.h>
#include "../../utils/math.h"

using namespace SAMRAI;

namespace cosmo
{
  
void scalar_ic_set_semianalytic_test(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  boost::shared_ptr<tbox::Database> cosmo_static_db);
}

#endif
