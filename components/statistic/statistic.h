#ifndef COSMO_STATISTIC_H
#define COSMO_STATISTIC_H

#include "../../cosmo_includes.h"
#include "../bssn/bssn.h"

using namespace SAMRAI;

namespace cosmo{
class CosmoStatistic
{
 public:

  CosmoStatistic(
    const tbox::Dimension& dim_in,
    boost::shared_ptr<tbox::Database> cosmo_statistic_db_in,
    std::ostream* l_stream_in);
  
  const tbox::Dimension& dim;
  boost::shared_ptr<tbox::Database> &cosmo_statistic_db;
  std::ostream* lstream;


  std::vector<std::string> conformal_avg_list;
  idx_t conformal_avg_interval;
  std::vector<idx_t> conformal_avg_idx;
  
  void output_conformal_avg(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    BSSN *bssn,
    idx_t weight_idx,
    idx_t step_num,
    real_t time);

  bool is_empty;
};

}
#endif
