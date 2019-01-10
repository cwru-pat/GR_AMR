/** @file bssn_ic.h
 * @brief Functions to set initial conditions for the BSSN class.
 * Functions are called via a config setting in VacuumSim
 * class.
 */

#ifndef COSMO_BSSN_ICS
#define COSMO_BSSN_ICS

#include "bssn.h"
#include "../../cosmo_includes.h"

using namespace SAMRAI;

#define INDEX(i,j,k) \
  (STENCIL_ORDER * (1 + (NX+ 2*STENCIL_ORDER) + (NX+ 2*STENCIL_ORDER) * (NY+ 2*STENCIL_ORDER))        \
   + (i) + (j) * (NX+ 2*STENCIL_ORDER) + (k) * (NX+ 2*STENCIL_ORDER) * (NY+ 2*STENCIL_ORDER))


namespace cosmo
{
  
void bssn_ic_awa_stability(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln, real_t A);
void bssn_ic_awa_linear_wave(BSSN * bssn);
void bssn_ic_awa_linear_wave(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln, real_t A, idx_t dir);
void bssn_ic_awa_linear_wave_desitter(BSSN * bssn);
void bssn_ic_awa_gauge_wave(BSSN * bssn);
void bssn_ic_awa_gauge_wave(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln, idx_t dir);
void bssn_ic_awa_shifted_gauge_wave(BSSN * bssn);
void bssn_ic_awa_shifted_gauge_wave(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln, idx_t dir);
void bssn_ic_static_blackhole(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln);

void bssn_ic_kerr_blackhole(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln);

void bssn_ic_ds_blackhole(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn);

void bssn_ic_kerr_BHL_CTT(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln,
  real_t M,
  real_t a,
  real_t K_c,
  real_t relaxation_tolerance,
  idx_t num_vcycles,
  idx_t max_depth,
  double l,
  double sigma);

void bssn_ic_static_BHL_CTT(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln,
  real_t M,
  real_t a,
  real_t K_c,
  real_t relaxation_tolerance,
  idx_t num_vcycles);

 
}

#endif
