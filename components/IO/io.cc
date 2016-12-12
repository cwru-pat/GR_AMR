#include "../../cosmo_includes.h"
#include "io.h"

using namespace SAMRAI

namespace cosmo{

CosmoIO:CosmoIO(
  const tbox::Dimension& dim_in,
  tbox::Database& cosmo_io_db_in
  std::ostream* l_stream_in = 0):
  dim(dim_in),
  cosmo_io_db(cosmo_io_db_in),
  lstream(l_stream_in),
{
  output_list = cosmo_sim_db->getStringVector("output_list");
  output_interval = cosmo_sim_db->getIntegerVector("output_interval");
  TBOX_ASSERT(output_list.size() == output_interval.size());
}

bool CosmoIO:packDerivedDataIntoDoubleBuffer(
  double* buffer,
  const hier::Patch& patch,
  const hier::Box& region,
  const std::string& variable_name,
  int depth_id,
  double simulation_time)
{
  NULL_USE(depth_id);
  NULL_USE(simulation_time);

   return true;
}

void CosmoIO::registerVariablesWithPlotter(
  appu::VisItDataWriter& visit_writer,
  idx_t step)
{
  
  for(idx_t i = 0; i < output_list.length(); i ++)
  {
    if(step % output_interval[i] != 0 ) continue;
    idx_t idx =
      variable_db->mapVariableAndContextToIndex(
        variable_db->getVariable(output_list[i]), variable_db->getContext("PREVIOUS"));
    
    if(idx < 0)
      TBOX_ERROR("Cannot find variable" <<*it<<" in output list!\n");
      
    visit_writer.registerPlotQuantity(
      output_list[i],
      "SCALAR",
      idx,
      0,
      1.0,
      "CELL");
  }
}

void CosmoIO::dumpData(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t step_num,
  real_t time)
{
  TBOX_ASSERT(visit_writer);
  
  visit_writer->writePlotData(
    patch_hierarchy,
    step_num,
    time);

  tbox::plog << "Wrote viz file for grid number "
             << step_num << '\n';
}


