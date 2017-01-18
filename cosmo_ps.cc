
#include "cosmo_includes.h"
#include "cosmo_ps.h"

/*
 * Headers for basic SAMRAI objects used in this code.
 */

using namespace SAMRAI;

namespace cosmo
{


CosmoPatchStrategy::CosmoPatchStrategy(
  const tbox::Dimension& dim_in,
  std::string object_name_in):
  xfer::RefinePatchStrategy(),
  dim(dim_in),
  object_name(object_name_in),
  is_time_dependent(false)
{

}


CosmoPatchStrategy::~CosmoPatchStrategy() {
}

void CosmoPatchStrategy::addTarget(idx_t idx)
{
  target_id_list.push_back(idx);
}    

}
