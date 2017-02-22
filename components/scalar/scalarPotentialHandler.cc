#include "scalarPotentialHandler.h"
#include "../../cosmo_types.h"
#include "../../utils/math.h"
#include <map>
#include <cmath>


using namespace SAMRAI;

namespace cosmo
{

/**
 * @brief Don't evolve anything
 * @return 0
 */
real_t scalarPotentialHandler::constant(
  BSSNData *bd, ScalarData *sd)
{
  return Lambda;
}

real_t scalarPotentialHandler::quadratic(
  BSSNData *bd, ScalarData *sd)
{
  return q_coef * pw2(sd->phi);
}
real_t scalarPotentialHandler::der_constant(
  BSSNData *bd, ScalarData *sd)
{
  return 0.0;
}

real_t scalarPotentialHandler::der_quadratic(
  BSSNData *bd, ScalarData *sd)
{
  return 2.0 * q_coef * sd->phi;
}

} // namespace cosmo


