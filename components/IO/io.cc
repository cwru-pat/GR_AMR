#include "../../cosmo_includes.h"
#include "io.h"

using namespace SAMRAI;

namespace cosmo{

CosmoIO::CosmoIO(
  const tbox::Dimension& dim_in,
  std::shared_ptr<tbox::Database> cosmo_io_db_in,
  std::ostream* l_stream_in):
  dim(dim_in),
  cosmo_io_db(cosmo_io_db_in),
  lstream(l_stream_in),
  is_empty(true)
{
  output_list = cosmo_io_db->getStringVector("output_list");
  output_interval = cosmo_io_db->getIntegerVector("output_interval");
#if !USE_COSMOTRACE
  output_null_geodesic = false;
  output_null_geodesic_interval = 0;
#endif
  TBOX_ASSERT(output_list.size() == output_interval.size());
}

/**
 * @brief  only need it when you need to output derived varibles(variables are not stored directly)
 */
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

  const int * lower = &region.lower()[0];
  const int * upper = &region.upper()[0];

  if(variable_name == "example")
  {
  
#pragma omp parallel for collapse(2)
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          int idx = (i-lower[0]) + (j - lower[1]) * (upper[0] - lower[0] + 1) +
            (k-lower[2]) * (upper[0] - lower[0] + 1) * (upper[1] - lower[1] + 1);
          double value = 0;

          buffer[idx] = idx;

        }
      }
    }
  }
  else
    TBOX_ERROR("The derived variable "<<variable_name<<" does not exist!\n");
  
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
        if(output_list[i] == "null_particles")
    {
#if !USE_COSMOTRACE
      TBOX_ERROR("Setting to output ray particle without turn on USE_COSMOTRACE!");
#else
      output_null_geodesic = true;
      output_null_geodesic_interval = output_interval[i];
#endif
      continue;
    }

    is_empty = 0;


    if(variable_db->getVariable(output_list[i]) == NULL)
    {
      tbox::plog<<"Cannot find variable" <<output_list[i]<<" in theoutput list!"
                <<" Registed as derived variable!\n";

      visit_writer.registerDerivedPlotQuantity(
        output_list[i],
        "SCALAR",
        (appu::VisDerivedDataStrategy *)this
      );
      
      continue;
    }
    
    idx_t idx =
      variable_db->mapVariableAndContextToIndex(
        variable_db->getVariable(output_list[i]), variable_db->getContext("ACTIVE"));    
      
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
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
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

#if USE_COSMOTRACE
void CosmoIO::dumpData(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  appu::VisItDataWriter& visit_writer,
  idx_t step_num,
  real_t time,
  Geodesic *ray,
  std::string hdf_filename_in)
{
  //TBOX_ASSERT(visit_writer);
  if(!is_empty)
  {
  
    visit_writer.writePlotData(
      hierarchy,
      step_num,
      time);

    tbox::plog << "Wrote viz file for grid number "
               << step_num << '\n';
  }

  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());

  if(step_num % output_null_geodesic_interval != 0)
    return;

  std::string dump_dirname = hdf_filename_in;
  dump_dirname += ".hdf/hdf_dump." + tbox::Utilities::intToString(step_num, 5);

  tbox::Utilities::recursiveMkdir(dump_dirname);
  
  std::string filename = dump_dirname + "/null_geodesic_proc_";
  filename += tbox::Utilities::intToString(mpi.getRank(), 5);
  filename += ".hdf";
  std::shared_ptr<tbox::HDFDatabase > hdf (new tbox::HDFDatabase("null_geodesic"));

  std::ifstream file(filename);

  // if(file)
  //   if(!remove(filename.c_str()))
  //     TBOX_ERROR("Failed to remove file "<<filename<<"\n");
  hdf->create(filename);
    
  if(!hdf->open(filename, true))
    TBOX_ERROR("Failed to open file "<<filename<<"\n");

  // int my_proc = mpi.getRank();
  // char temp_buf[20];
  // sprintf(temp_buf, "processor.%05d", my_proc);
  // std::shared_ptr<tbox::Database> processor_HDFGroup(
  //   hdf->putDatabase(std::string(temp_buf)));

  
  int particle_cnt = 0;

  std::vector<double> all_info;
  
  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
  {
    std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));
    // put all particles to correct cell after advancing it
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end(); ++ip)
    {
      std::shared_ptr<hier::Patch> patch(*ip);
      ray->initPData(patch);
      const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      const double *dx = &patch_geom->getDx()[0];

      const hier::Box& box = patch->getBox();

      pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iter(*ray->pc_pdata, true);
      pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iterend(*ray->pc_pdata, false);

      for(;iter != iterend; iter++)
      {
        ParticleContainer &id = *iter;

        if(!box.contains(id.idx))
          continue;
        
        //        double p_info[PARTICLE_NUMBER_OF_STATES + PARTICLE_INT_PROPERTIES];
        
        for(std::list<RKParticle>::iterator it=id.p_list.begin();
            it != id.p_list.end(); it++)
        {
          for(int i = 0; i < PARTICLE_NUMBER_OF_STATES; i ++)
            all_info.push_back((*it).x_a[i]);
          for(int i = 0; i < PARTICLE_REAL_PROPERTIES; i ++)
            all_info.push_back((*it).rp[i]);                    
          for(int i = 0; i < PARTICLE_INT_PROPERTIES; i ++)
            all_info.push_back((*it).ip[i]);
          
          particle_cnt++;
        }
      }
    }
  }
  for(int i = 0; i < all_info.size(); i++)
    std::cout<<all_info[i]<<" ";
  std::cout<<"\n";

  if(all_info.size() > 0)
  {
    hdf->putDoubleArray("null_geodesic", &all_info[0], all_info.size());
   
    hdf->putInteger("num_of_particels", particle_cnt);
  }
  hdf->close();
}
#endif
  
/**
 * @brief print patch for debuging
 */
void CosmoIO::printPatch(
  const std::shared_ptr<hier::Patch> & patch,
  std::ostream &os,
  idx_t idx)
{
  std::shared_ptr<pdat::CellData<double> > pd(
    SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
      patch->getPatchData(idx)));
  const hier::Box& ghost_box = pd->getGhostBox();
  pd->print(ghost_box, 0, os, 9);
}
  
}
