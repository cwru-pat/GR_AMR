#ifndef COSMO_IO_H
#define COSMO_IO_H

#include "../../cosmo_includes.h"

using namespace SAMRAI

namespace cosmo{
class CosmoIO:public appu::VisDerivedDataStrategy
{
 public:

  CosmoIO(
    const tbox::Dimension& dim_in,
    tbox::Database& cosmo_io_db_in
    std::ostream* l_stream_in = 0);
  
  const tbox::Dimension& dim_in;
  boost::shared_ptr<tbox::Database> cosmo_io_db;
  std::ostream* lstream;


  std::vector<std::string> output_list;
  std::vector<std::idx_t> output_interval;
  
  virtual bool
   packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& variable_name,
      int depth_id,
      double simulation_time) const;

  void registerVariablesWithPlotter(
    appu::VisItDataWriter& visit_writer,
    idx_t step);

  void dumpData(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    idx_t step_num,
    real_t time);

  
};

}
#endif
