#ifndef COSMO_IO_H
#define COSMO_IO_H

#include "../../cosmo_includes.h"

using namespace SAMRAI;

namespace cosmo{
class CosmoIO:public appu::VisDerivedDataStrategy
{
 public:

  CosmoIO(
    const tbox::Dimension& dim_in,
    std::shared_ptr<tbox::Database> cosmo_io_db_in,
    std::ostream* l_stream_in);
  
  const tbox::Dimension& dim;
  std::shared_ptr<tbox::Database> &cosmo_io_db;
  std::ostream* lstream;


  std::vector<std::string> output_list;
  std::vector<idx_t> output_interval;
  
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
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    appu::VisItDataWriter& visit_writer,
    idx_t step_num,
    real_t time);
  void printPatch(
    const std::shared_ptr<hier::Patch> & patch,
    std::ostream &os,
    idx_t idx);

  bool is_empty;
};

}
#endif
