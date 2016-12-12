#ifndef COSMO_BD_H
#define COSMO_BD_H

#include "cosmo_includes.h"
#include "cosmo_ps.h"

/*
 * Headers for basic SAMRAI objects used in this code.
 */

using namespace SAMRAI;

namespace cosmo
{

CosmoBD: CosmoPatchStrategy()
{
  
}

CosmoBD: ~CosmoPatchStrategy()
{
  
}

hier::IntVector CosmoPatchStrategy:getRefineOpStencilWidth(
      const tbox::Dimension& dim)
{
  return hier::IntVector::getZero(dim);
}
void CosmoPatchStrategy:addTarget(idx_t idx)
{
  target_id_list.push_back(idx);
}    

}
