/** @file static_ic.h
 * @brief Functions to set initial conditions for the Static (dust) class.
 * Functions should be made callable via a config setting in the DustSim
 * class.
 */

#ifndef COSMO_STATIC_ICS
#define COSMO_STATIC_ICS

#include "../../cosmo_includes.h"
#include "static.h"
#include "../../ICs/ICs.h"
#include "../../ICs/ICs_data.h"
#include "../../utils/math.h"

using namespace SAMRAI;

namespace cosmo
{
  
void dust_ic_set_random(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, std::shared_ptr<tbox::Database> cosmo_static_db);

}

#endif
