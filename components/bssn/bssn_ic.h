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

namespace cosmo
{
  
void bssn_ic_awa_stability(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln, real_t A);
void bssn_ic_awa_linear_wave(BSSN * bssn);
void bssn_ic_awa_linear_wave(BSSN * bssn, real_t A, int dir);
void bssn_ic_awa_linear_wave_desitter(BSSN * bssn);
void bssn_ic_awa_gauge_wave(BSSN * bssn);
void bssn_ic_awa_gauge_wave(BSSN * bssn, int dir);
void bssn_ic_awa_shifted_gauge_wave(BSSN * bssn);
void bssn_ic_awa_shifted_gauge_wave(BSSN * bssn, int dir);
void bssn_ic_static_blackhole(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln);
 
}

#endif
