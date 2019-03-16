#include "../../cosmo_includes.h"
#include "statistic.h"

using namespace SAMRAI;

namespace cosmo{

CosmoStatistic::CosmoStatistic(
  const tbox::Dimension& dim_in,
  std::shared_ptr<tbox::Database> cosmo_statistic_db_in,
  std::ostream* l_stream_in):
  dim(dim_in),
  cosmo_statistic_db(cosmo_statistic_db_in),
  lstream(l_stream_in),
  is_empty(true)
{
  conformal_avg_interval = cosmo_statistic_db->getIntegerWithDefault("conformal_avg_interval",0);
  expansion_info_interval = cosmo_statistic_db->getIntegerWithDefault("expansion_info_interval",0);

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

real_t CosmoStatistic::calculate_conformal_avg(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  BSSN *bssn,
  idx_t weight_idx,
  idx_t field_idx,
  //calculating conformal avg only on boundary
  bool only_on_bd)
{
  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;
  const double * dx = &grid_geometry.getDx()[0];

  real_t L[3];
  const double * domain_lower = &grid_geometry.getXLower()[0];
  const double * domain_upper = &grid_geometry.getXUpper()[0];

  for(int i = 0 ; i < 3; i++)
    L[i] = domain_upper[i] - domain_lower[i];

  int output_num = static_cast<idx_t>(conformal_avg_list.size());
  double conformal_avg = 0;

  arr_t * arrays = new arr_t[output_num];

  double tot_vol = 0;

  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln ++)
  {
    std::shared_ptr <hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

    
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
    {
      const std::shared_ptr<hier::Patch> & patch = *pit;

      const hier::Box& box = patch->getBox();

      const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      bssn->initPData(patch);

      bssn->initMDA(patch);


      std::shared_ptr<pdat::CellData<double> > weight(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_idx)));
      
      arr_t weight_array =
        pdat::ArrayDataAccess::access<DIM, double>(
          weight->getArrayData());


      std::shared_ptr<pdat::CellData<double> > var_cell_data(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(field_idx)));

      arr_t array = pdat::ArrayDataAccess::access<DIM, double>(
        var_cell_data->getArrayData());
      
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];
      
      const double *dx = &patch_geom->getDx()[0];

      
#pragma omp parallel for collapse(2) reduction(+:tot_vol, conformal_avg)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            if(weight_array(i,j,k) > 0)
            {
              tot_vol += weight_array(i, j, k) * 1.0 / pw3( + 1.0);
              conformal_avg +=
                1.0 / pw3(+ 1.0) * array(i, j, k)
                * weight_array(i, j, k);
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
    mpi.AllReduce(&conformal_avg, 1, MPI_SUM);
  }

  return conformal_avg / tot_vol;
  
}

void CosmoStatistic::output_expansion_info(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    BSSN *bssn,
    idx_t weight_idx,
    idx_t step_num,
    real_t time,
    real_t min_radius)
{
  
}
  
void CosmoStatistic::output_conformal_avg(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    BSSN *bssn,
    idx_t weight_idx,
    idx_t step_num,
    real_t time)
{
  if(conformal_avg_interval == 0 || step_num % conformal_avg_interval != 0)
    return;
  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;
  const double * dx = &grid_geometry.getDx()[0];

  int output_num = static_cast<idx_t>(conformal_avg_list.size());
  double *conformal_avg = new double[output_num];

  arr_t * arrays = new arr_t[output_num];

  double tot_vol = 0;

  double max_r = 0;

  int r_idx = -1;

  for(int i = 0; i < output_num; i++)
    if(conformal_avg_list[i] == "DIFFr")
    {
      r_idx = i;
      break;
    }
  if(r_idx == -1)
    TBOX_ERROR("No DIFFr as statistics input, cannot output max_rho / <rho>");
  
  for(int i = 0; i < output_num; i++) conformal_avg[i] = 0;
  
  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln ++)
  {
    std::shared_ptr <hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

    
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
    {
      const std::shared_ptr<hier::Patch> & patch = *pit;

      const hier::Box& box = patch->getBox();

      const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      bssn->initPData(patch);

      bssn->initMDA(patch);


      std::shared_ptr<pdat::CellData<double> > weight(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_idx)));
      
      arr_t weight_array =
        pdat::ArrayDataAccess::access<DIM, double>(
          weight->getArrayData());


      for(int i = 0; i < output_num; i ++)
      {
        std::shared_ptr<pdat::CellData<double> > var_cell_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch->getPatchData(conformal_avg_idx[i])));

        arrays[i] = pdat::ArrayDataAccess::access<DIM, double>(
          var_cell_data->getArrayData());
      }
      
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];
      
      const double *dx = &patch_geom->getDx()[0];

      

#pragma omp parallel for collapse(2) reduction(+:tot_vol) reduction(max:max_r)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            
            if(weight_array(i,j,k) > 0)
            {
              tot_vol += weight_array(i, j, k) * 1.0 / pw3( + 1.0);

              max_r = std::max(max_r, arrays[r_idx](i, j, k));

              #pragma omp critical
              {
                for(int var_id = 0; var_id < output_num; var_id++)
                {
                  conformal_avg[var_id] +=
                    1.0 / pw3( + 1.0) * arrays[var_id](i, j, k)
                    * weight_array(i, j, k);
                }
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
    mpi.AllReduce(&max_r, 1, MPI_MAX);
    mpi.AllReduce(&conformal_avg[0], output_num, MPI_SUM);
  }

  (*lstream)<<"Total vol is " << tot_vol<<"\n";
  
  for(int var_id = 0; var_id < output_num; var_id++)
  {
    (*lstream)<<"Conformal Avg for "<<conformal_avg_list[var_id]
            <<" is "
            <<std::setprecision(12)<<conformal_avg[var_id] / tot_vol<<"\n";
  }

  (*lstream)<<"rho_max / rho_avg is "<<max_r / (conformal_avg[r_idx]/tot_vol)<<"\n";
  
}

// only do statistics on the sphere larger than min_radius
void CosmoStatistic::output_conformal_avg(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    BSSN *bssn,
    idx_t weight_idx,
    idx_t step_num,
    real_t time,
    real_t min_radius)
{
  if(conformal_avg_interval == 0 || step_num % conformal_avg_interval != 0)
    return;

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  const double * domain_lower = &grid_geometry.getXLower()[0];
  const double * domain_upper = &grid_geometry.getXUpper()[0];

  double L[DIM];
  for(int i = 0 ; i < DIM; i++)
    L[i] = domain_upper[i] - domain_lower[i];

  const double * dx = &grid_geometry.getDx()[0];

  int output_num = static_cast<idx_t>(conformal_avg_list.size());
  double *conformal_avg = new double[output_num];

  arr_t * arrays = new arr_t[output_num];

  double tot_vol = 0;

  for(int i = 0; i < output_num; i++) conformal_avg[i] = 0;
  
  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln ++)
  {
    std::shared_ptr <hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

    
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
    {
      const std::shared_ptr<hier::Patch> & patch = *pit;

      const hier::Box& box = patch->getBox();

      const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      bssn->initPData(patch);

      bssn->initMDA(patch);


      std::shared_ptr<pdat::CellData<double> > weight(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_idx)));
      
      arr_t weight_array =
        pdat::ArrayDataAccess::access<DIM, double>(
          weight->getArrayData());


      for(int i = 0; i < output_num; i ++)
      {
        std::shared_ptr<pdat::CellData<double> > var_cell_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
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

            double xx = (dx[0] * ((real_t)i + 0.5)) - L[0] / 2.0 ;
            double yy = (dx[1] * ((real_t)j + 0.5)) - L[1] / 2.0 ;
            double zz = (dx[2] * ((real_t)k + 0.5)) - L[2] / 2.0 ;
            double r = sqrt(pw2(xx) + pw2(yy) + pw2(zz));

            if(weight_array(i,j,k) > 0 && r > min_radius)
            {
              tot_vol += weight_array(i, j, k) * 1.0 / pw3( + 1.0);
              #pragma omp critical
              {
                for(int var_id = 0; var_id < output_num; var_id++)
                  conformal_avg[var_id] +=
                    1.0 / pw3( + 1.0) * arrays[var_id](i, j, k)
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
