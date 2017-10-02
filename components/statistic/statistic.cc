#include "../../cosmo_includes.h"
#include "statistic.h"

using namespace SAMRAI;

namespace cosmo{

CosmoStatistic::CosmoStatistic(
  const tbox::Dimension& dim_in,
  boost::shared_ptr<tbox::Database> cosmo_statistic_db_in,
  std::ostream* l_stream_in):
  dim(dim_in),
  cosmo_statistic_db(cosmo_statistic_db_in),
  lstream(l_stream_in),
  is_empty(true)
{
  conformal_avg_interval = cosmo_statistic_db->getIntegerWithDefault("conformal_avg_interval",0);

  if(conformal_avg_interval > 0)
    conformal_avg_list = cosmo_statistic_db->getStringVector("conformal_avg_list");

  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

  for(idx_t i = 0; i < static_cast<idx_t>(conformal_avg_list.size()); i ++)
  {
    idx_t idx =
      variable_db->mapVariableAndContextToIndex(
        variable_db->getVariable(conformal_avg_list[i]), variable_db->getContext("ACTIVE"));
    
    if(idx < 0)
      TBOX_ERROR("Cannot find variable" <<conformal_avg_list[i]<<" in output list!\n");
    
    conformal_avg_idx.push_back(idx);
    
  }
  
}

void CosmoStatistic::output_conformal_avg(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    BSSN *bssn,
    idx_t weight_idx,
    idx_t step_num,
    real_t time)
{
  if(conformal_avg_interval == 0 || step_num % conformal_avg_interval != 0)
    return;
  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;
  const double * dx = &grid_geometry.getDx()[0];

  int output_num = static_cast<idx_t>(conformal_avg_list.size());
  double *conformal_avg = new double[output_num];

  arr_t * arrays = new arr_t[output_num];

  double tot_vol = 0;

  for(int i = 0; i < output_num; i++) conformal_avg[i] = 0;
  
  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln ++)
  {
    boost::shared_ptr <hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

    
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
    {
      const boost::shared_ptr<hier::Patch> & patch = *pit;

      const hier::Box& box = patch->getBox();

      const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      bssn->initPData(patch);

      bssn->initMDA(patch);


      boost::shared_ptr<pdat::CellData<double> > weight(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_idx)));
      
      arr_t weight_array =
        pdat::ArrayDataAccess::access<DIM, double>(
          weight->getArrayData());

      arr_t DIFFchi_array = bssn->DIFFchi_a;

      for(int i = 0; i < output_num; i ++)
      {
        boost::shared_ptr<pdat::CellData<double> > var_cell_data(
          BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            patch->getPatchData(conformal_avg_idx[i])));

        arrays[i] = pdat::ArrayDataAccess::access<DIM, double>(
          var_cell_data->getArrayData());
      }
      
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];
      
      const double *dx = &patch_geom->getDx()[0];

      

#pragma omp parallel for collapse(2) reduction(+:tot_vol)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            
            if(weight_array(i,j,k) > 0)
            {
              tot_vol += weight_array(i, j, k) * 1.0 / pw3(DIFFchi_array(i, j, k) + 1.0);
              #pragma omp critical
              {
                for(int var_id = 0; var_id < output_num; var_id++)
                  conformal_avg[var_id] +=
                    1.0 / pw3(DIFFchi_array(i, j, k) + 1.0) * arrays[var_id](i, j, k)
                    * weight_array(i, j, k);
              }
            }
          }
        }
      }
    }
  }
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  mpi.Barrier();
  if (mpi.getSize() > 1) {
    mpi.AllReduce(&tot_vol, 1, MPI_SUM);
    mpi.AllReduce(&conformal_avg[0], output_num, MPI_SUM);
  }

  (*lstream)<<"Total vol is " << tot_vol<<"\n";
  
  for(int var_id = 0; var_id < output_num; var_id++)
  {
    (*lstream)<<"Conformal Avg for "<<conformal_avg_list[var_id]
            <<" is "
            <<std::setprecision(12)<<conformal_avg[var_id] / tot_vol<<"\n";
  }

}


}
