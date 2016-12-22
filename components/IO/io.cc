#include "../../cosmo_includes.h"
#include "io.h"

using namespace SAMRAI;

namespace cosmo{

CosmoIO::CosmoIO(
  const tbox::Dimension& dim_in,
  boost::shared_ptr<tbox::Database> cosmo_io_db_in,
  std::ostream* l_stream_in):
  dim(dim_in),
  cosmo_io_db(cosmo_io_db_in),
  lstream(l_stream_in),
  is_empty(true)
{
  output_list = cosmo_io_db->getStringVector("output_list");
  output_interval = cosmo_io_db->getIntegerVector("output_interval");
  TBOX_ASSERT(output_list.size() == output_interval.size());
}

bool CosmoIO::packDerivedDataIntoDoubleBuffer(
  double* buffer,
  const hier::Patch& patch,
  const hier::Box& region,
  const std::string& variable_name,
  int depth_id,
  double simulation_time) const
{
  NULL_USE(depth_id);
  NULL_USE(simulation_time);

   return true;
}

void CosmoIO::registerVariablesWithPlotter(
  appu::VisItDataWriter& visit_writer,
  idx_t step)
{
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  is_empty = 1;
  for(idx_t i = 0; i < static_cast<idx_t>(output_list.size()); i ++)
  {
    if(step % output_interval[i] != 0 ) continue;
    is_empty = 0;
    idx_t idx =
      variable_db->mapVariableAndContextToIndex(
        variable_db->getVariable(output_list[i]), variable_db->getContext("ACTIVE"));
    
    if(idx < 0)
      TBOX_ERROR("Cannot find variable" <<output_list[i]<<" in output list!\n");
      
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
  appu::VisItDataWriter& visit_writer,
  idx_t step_num,
  real_t time)
{
  //TBOX_ASSERT(visit_writer);

  if(is_empty) return;
  
  visit_writer.writePlotData(
    hierarchy,
    step_num,
    time);

  tbox::plog << "Wrote viz file for grid number "
             << step_num << '\n';
}

void CosmoIO::printPatch(
  const boost::shared_ptr<hier::Patch> & patch,
  std::ostream &os,
  idx_t idx)
{
  boost::shared_ptr<pdat::CellData<double> > pd(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch->getPatchData(idx)));
  const hier::Box& ghost_box = pd->getGhostBox();
  pd->print(ghost_box, 0, os, PRINT_PRECISION);
}
  
}
