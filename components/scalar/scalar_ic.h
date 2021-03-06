#ifndef COSMO_SCALAR_ICS
#define COSMO_SCALAR_ICS

#include "../../cosmo_includes.h"
#include "scalar.h"
#include <fftw3.h>
#include <zlib.h>
#include "../../utils/math.h"
#include "../bssn/bssn_data.h"
#include "scalar_data.h"
#include <fstream>

using namespace SAMRAI;

#define INDEX(i,j,k) \
  (STENCIL_ORDER * (1 + (NX+ 2*STENCIL_ORDER) + (NX+ 2*STENCIL_ORDER) * (NY+ 2*STENCIL_ORDER))        \
   + (i) + (j) * (NX+ 2*STENCIL_ORDER) + (k) * (NX+ 2*STENCIL_ORDER) * (NY+ 2*STENCIL_ORDER))

#define LOOP3()  \
  for(int i=0; i<NX; ++i) \
    for(int j=0; j<NY; ++j) \
      for(int k=0; k<NZ; ++k)
#define NP_INDEX(i,j,k) ((NZ)*(NY)*(i) + (NZ)*(j) + (k))
#define FFT_NP_INDEX(i,j,k) ((NZ/2+1)*NY*(i) + (NZ/2+1)*(j) + (k))


namespace cosmo
{
  
void scalar_ic_set_semianalytic_test(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  std::shared_ptr<tbox::Database> cosmo_static_db);

 bool scalar_ic_set_scalar_collapse(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  std::shared_ptr<tbox::Database> cosmo_scalar_db);

 bool scalar_ic_set_oscillon(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  std::shared_ptr<tbox::Database> cosmo_scalar_db);
 
 bool scalar_ic_set_scalar_gaussian_collapse(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  std::shared_ptr<tbox::Database> cosmo_scalar_db);

 
 bool scalar_ic_set_scalar_collapse_sommerfield(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  std::shared_ptr<tbox::Database> cosmo_scalar_db);
 bool scalar_ic_set_periodic_fast_collapse_test(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   idx_t ln, BSSN * bssn, Scalar * scalar,
   std::shared_ptr<tbox::Database> cosmo_scalar_db);
bool scalar_ic_set_scalar_gaussian_random(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  std::shared_ptr<tbox::Database> cosmo_scalar_db);

 
 inline bool exist(const std::string& name)
 {
   std::ifstream file(name);
   if(!file)            // If the file was not found, then file is 0, i.e. !file=1 or true.
     return false;    // The file was not found.
   else                 // If the file was found, then file is non-0.
     return true;     // The file was found.
 }
 
}

#endif
