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
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  BSSN *bssn,
  idx_t weight_idx,
  idx_t field_idx,
  //calculating conformal avg only on boundary
  bool only_on_bd)
{
  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
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

      arr_t gamma11 = bssn->DIFFgamma11_a;
      arr_t gamma12 = bssn->DIFFgamma12_a;
      arr_t gamma13 = bssn->DIFFgamma13_a;
      arr_t gamma22 = bssn->DIFFgamma22_a;
      arr_t gamma23 = bssn->DIFFgamma23_a;
      arr_t gamma33 = bssn->DIFFgamma33_a;

      
      boost::shared_ptr<pdat::CellData<double> > var_cell_data(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
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
            if(only_on_bd)
            {
              if(i == 0 && j == 0)
              {
                tot_vol +=  dx[2] * (gamma33(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0);
                conformal_avg += dx[2] * (gamma33(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0)
                  * array(i, j, k);
              }
              else if(i == 0 && j == (round(L[1]/dx[1])) - 1)
              {
                tot_vol += dx[2] * (gamma33(i, j, k) + 1.0)  / (DIFFchi_array(i, j, k) + 1.0);
                conformal_avg += dx[2] * (gamma33(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0)
                  * array(i, j, k);
              }
              else if(j == 0 && i == round(L[0]/dx[0]) - 1)
              {
                tot_vol += dx[2] * (gamma33(i, j, k) + 1.0)  / (DIFFchi_array(i, j, k) + 1.0);
                conformal_avg += dx[2] * (gamma33(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0)
                  * array(i, j, k);

              }
              else if(i == round(L[0]/dx[0]) - 1 && j == (round(L[1]/dx[1])) - 1)
              {
                tot_vol += dx[2] * (gamma33(i, j, k) + 1.0)  / (DIFFchi_array(i, j, k) + 1.0);
                conformal_avg += dx[2] * (gamma33(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0)
                  * array(i, j, k);
              }
              if(i == 0 && k == 0)
              {
                tot_vol +=  dx[1] * (gamma22(i, j, k) + 1.0)  / (DIFFchi_array(i, j, k) + 1.0);
                conformal_avg += dx[1] * (gamma22(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0)
                  * array(i, j, k);
              }
              else if(i == 0 && k == (round(L[2]/dx[2])) - 1)
              {
                tot_vol += dx[1] * (gamma22(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0);
                conformal_avg += dx[1] * (gamma22(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0)
                  * array(i, j, k);
              }
              else if(k == 0 && i == round(L[0]/dx[0]) - 1)
              {
                tot_vol += dx[1] * (gamma22(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0);
                conformal_avg += dx[1] * (gamma22(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0)
                  * array(i, j, k);
              }
              else if(i == round(L[0]/dx[0]) - 1 && k == (round(L[2]/dx[2])) - 1)
              {
                tot_vol += dx[1] * (gamma22(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0);
                conformal_avg += dx[1] * (gamma22(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0)
                  * array(i, j, k);
              }
                            
              if(k == 0 && j == 0)
              {
                tot_vol +=  dx[0] * (gamma11(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0);
                conformal_avg += dx[0] * (gamma11(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0)
                  * array(i, j, k);
              }
              else if(j == 0 && k == (round(L[2]/dx[2])) - 1)
              {
                tot_vol += dx[0] * (gamma11(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0);
                conformal_avg += dx[0] * (gamma11(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0)
                  * array(i, j, k);
              }
              else if(k == 0 && j == (round(L[1]/dx[1])) - 1)
              {
                tot_vol += dx[0] * (gamma11(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0);
                conformal_avg += dx[0] * (gamma11(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0)
                  * array(i, j, k);

              }
              else if(j == (round(L[1]/dx[1])) - 1 && k == (round(L[2]/dx[2])) - 1)
              {
                tot_vol += dx[0] * (gamma11(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0);
                conformal_avg += dx[0] * (gamma11(i, j, k) + 1.0) / (DIFFchi_array(i, j, k) + 1.0)
                  * array(i, j, k);
              }

              continue;
            }
            if(weight_array(i,j,k) > 0)
            {
              tot_vol += weight_array(i, j, k) * 1.0 / pw3(DIFFchi_array(i, j, k) + 1.0);
              conformal_avg +=
                1.0 / pw3(DIFFchi_array(i, j, k) + 1.0) * array(i, j, k)
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
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    BSSN *bssn,
    idx_t weight_idx,
    idx_t step_num,
    real_t time)
{
  if( expansion_info_interval == 0 || step_num % expansion_info_interval != 0)
    return;
  
  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;
  const double * dx = &grid_geometry.getDx()[0];

  const double * domain_lower = &grid_geometry.getXLower()[0];
  const double * domain_upper = &grid_geometry.getXUpper()[0];

  int coarsest_boxes[3];

  for(int i = 0; i < 3; i++)
    coarsest_boxes[i] = round((domain_upper[i] - domain_lower[i]) / dx[i] );
  
  double tot_vol = 0;
  double face_area[6] = {0}; //(0, y, z), (L, y, z), (x, 0, z), (x, L, z), (x, y, 0), (x, y, L)
  double K_face[6] = {0};
  
  //(0, 0, z), (0, L, z), (L, 0, z), (L, L, z)
  //(0, y, 0), (0, y, L), (L, 0, z), (L, y, L)
  //(x, 0, 0), (x, 0, L), (x, L, 0), (x, L, L)
  double edge_length[12] = {0};
  double K_edge[12] = {0};

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

      arr_t DIFFchi = bssn->DIFFchi_a;
      arr_t K = bssn->DIFFK_a;
      arr_t gamma11 = bssn->DIFFgamma11_a;
      arr_t gamma12 = bssn->DIFFgamma12_a;
      arr_t gamma13 = bssn->DIFFgamma13_a;
      arr_t gamma22 = bssn->DIFFgamma22_a;
      arr_t gamma23 = bssn->DIFFgamma23_a;
      arr_t gamma33 = bssn->DIFFgamma33_a;

      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];
      
      const double *dx = &patch_geom->getDx()[0];


      //#pragma omp parallel for collapse(2) reduction(+:tot_vol)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {           
            if(weight_array(i,j,k) > 0)
            {
              tot_vol += weight_array(i, j, k) * 1.0 / pw3(DIFFchi(i, j, k) + 1.0);

              if(i == 0)
              {
                double sq_det2 = ((gamma22(i, j, k) + 1.0) * (gamma33(i, j, k) + 1.0)
                                  - pw2(gamma23(i, j, k))) / pw2(DIFFchi(i, j, k) + 1.0);
                face_area[0] += dx[1] * dx[2] * sq_det2;
                K_face[0] += dx[1] * dx[2] * sq_det2 * K(i, j, k);
              }
              else if(i == (coarsest_boxes[0] << ln) - 1)
              {
                double sq_det2 = ( (gamma22(i, j, k) + 1.0)* (gamma33(i, j, k)  + 1.0)
                                   - pw2(gamma23(i, j, k))) / pw2(DIFFchi(i, j, k) + 1.0);
                face_area[1] += dx[1] * dx[2] * sq_det2;
                K_face[1] += dx[1] * dx[2] * sq_det2 * K(i, j, k);
              }
              if(j == 0)
              {
                double sq_det2 = ( (gamma11(i, j, k) + 1.0)* (gamma33(i, j, k) + 1.0)
                                   - pw2(gamma13(i, j, k))) / pw2(DIFFchi(i, j, k) + 1.0);
                face_area[2] += dx[0] * dx[2] * sq_det2;
                K_face[2] += dx[0] * dx[2] * sq_det2 * K(i, j, k);
              }
              else if(j == (coarsest_boxes[1] << ln) - 1)
              {
                double sq_det2 = ((gamma11(i, j, k) + 1.0) * (gamma33(i, j, k) + 1.0)
                                  - pw2(gamma13(i, j, k))) / pw2(DIFFchi(i, j, k) + 1.0);
                face_area[3] += dx[0] * dx[2] * sq_det2;
                K_face[3] += dx[0] * dx[2] * sq_det2 * K(i, j, k);
              }
              if(k == 0)
              {
                double sq_det2 = ((gamma11(i, j, k) + 1.0) * (gamma22(i, j, k) + 1.0)
                                  - pw2(gamma12(i, j, k))) / pw2(DIFFchi(i, j, k) + 1.0);
                face_area[4] += dx[0] * dx[1] * sq_det2;
                K_face[4] += dx[1] * dx[0] * sq_det2 * K(i, j, k);
              }
              else if(k == (coarsest_boxes[2] << ln) - 1)
              {
                double sq_det2 = ((gamma11(i, j, k) + 1.0) * (gamma22(i, j, k) + 1.0)
                                  - pw2(gamma12(i, j, k))) / pw2(DIFFchi(i, j, k) + 1.0);
                face_area[5] += dx[0] * dx[1] * sq_det2;
                K_face[5] += dx[1] * dx[0] * sq_det2 * K(i, j, k);
              }

              if(i == 0 && j == 0)
              {
                edge_length[0] +=  dx[2] * (gamma33(i, j, k) + 1.0) / (DIFFchi(i, j, k) + 1.0);
                K_edge[0] += K(i, j, k) * dx[2] * (gamma33(i, j, k) + 1.0) / (DIFFchi(i, j, k) + 1.0);
              }
              else if(i == 0 && j == (coarsest_boxes[1] << ln) - 1)
              {
                edge_length[1] += dx[2] * (gamma33(i, j, k) + 1.0)  / (DIFFchi(i, j, k) + 1.0);
                K_edge[1] += K(i, j, k) * dx[2] * (gamma33(i, j, k) + 1.0) / (DIFFchi(i, j, k) + 1.0);
              }
              else if(j == 0 && i == (coarsest_boxes[0] << ln) - 1)
              {
                edge_length[2] += dx[2] * (gamma33(i, j, k) + 1.0)  / (DIFFchi(i, j, k) + 1.0);
                K_edge[2] += K(i, j, k) * dx[2] * (gamma33(i, j, k) + 1.0) / (DIFFchi(i, j, k) + 1.0);
              }
              else if(i == (coarsest_boxes[0] << ln) - 1 && j == (coarsest_boxes[1] << ln) - 1)
              {
                edge_length[3] += dx[2] * (gamma33(i, j, k) + 1.0)  / (DIFFchi(i, j, k) + 1.0);
                K_edge[3] += K(i, j, k) * dx[2] * (gamma33(i, j, k) + 1.0) / (DIFFchi(i, j, k) + 1.0);
              }

              if(i == 0 && k == 0)
              {
                edge_length[4] +=  dx[1] * (gamma22(i, j, k) + 1.0)  / (DIFFchi(i, j, k) + 1.0);
                K_edge[4] += K(i, j, k) * dx[1] * (gamma22(i, j, k) + 1.0)  / (DIFFchi(i, j, k) + 1.0);
              }
              else if(i == 0 && k == (coarsest_boxes[2] << ln) - 1)
              {
                edge_length[5] += dx[1] * (gamma22(i, j, k) + 1.0) / (DIFFchi(i, j, k) + 1.0);
                K_edge[5] += K(i, j, k) * dx[1] * (gamma22(i, j, k) + 1.0)  / (DIFFchi(i, j, k) + 1.0);
              }
              else if(k == 0 && i == (coarsest_boxes[0] << ln) - 1)
              {
                edge_length[6] += dx[1] * (gamma22(i, j, k) + 1.0) / (DIFFchi(i, j, k) + 1.0);
                K_edge[6] += K(i, j, k) * dx[1] * (gamma22(i, j, k) + 1.0)  / (DIFFchi(i, j, k) + 1.0);
              }
              else if(i == (coarsest_boxes[0] << ln) - 1 && k == (coarsest_boxes[2] << ln) - 1)
              {
                edge_length[7] += dx[1] * (gamma22(i, j, k) + 1.0) / (DIFFchi(i, j, k) + 1.0);
                K_edge[7] += K(i, j, k) * dx[1] * (gamma22(i, j, k) + 1.0)  / (DIFFchi(i, j, k) + 1.0);
              }
                            
              if(k == 0 && j == 0)
              {
                edge_length[8] +=  dx[0] * (gamma11(i, j, k) + 1.0) / (DIFFchi(i, j, k) + 1.0);
                K_edge[8] += K(i, j, k) * dx[0] * (gamma11(i, j, k) + 1.0)  / (DIFFchi(i, j, k) + 1.0);
              }
              else if(j == 0 && k == (coarsest_boxes[2] << ln) - 1)
              {
                edge_length[9] += dx[0] * (gamma11(i, j, k) + 1.0) / (DIFFchi(i, j, k) + 1.0);
                K_edge[9] += K(i, j, k) * dx[0] * (gamma11(i, j, k) + 1.0)  / (DIFFchi(i, j, k) + 1.0);
              }
              else if(k == 0 && j == (coarsest_boxes[1] << ln) - 1)
              {
                edge_length[10] += dx[0] * (gamma11(i, j, k) + 1.0) / (DIFFchi(i, j, k) + 1.0);
                K_edge[10] += K(i, j, k) * dx[0] * (gamma11(i, j, k) + 1.0)  / (DIFFchi(i, j, k) + 1.0);
              }
              else if(j == (coarsest_boxes[1] << ln) - 1 && k == (coarsest_boxes[2] << ln) - 1)
              {
                edge_length[11] += dx[0] * (gamma11(i, j, k) + 1.0) / (DIFFchi(i, j, k) + 1.0);
                K_edge[11] += K(i, j, k) * dx[0] * (gamma11(i, j, k) + 1.0)  / (DIFFchi(i, j, k) + 1.0);
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
    mpi.AllReduce(&face_area[0], 6, MPI_SUM);
    mpi.AllReduce(&edge_length[0], 12, MPI_SUM);
    mpi.AllReduce(&K_face[0], 6, MPI_SUM);
    mpi.AllReduce(&K_edge[0], 12, MPI_SUM);

  }

  (*lstream).precision(7);
  (*lstream)<<"**********Start outputing expansion info****** " <<"\n";

  (*lstream)<<"Total vol is " <<tot_vol<<"\n";

  (*lstream)<<"6 faces areas are ";
  for(int i = 0; i < 6; i++)
    (*lstream)<<face_area[i]<<" ";
  (*lstream)<<"\n";

  (*lstream)<<"K on 6 faces areas are ";
  for(int i = 0; i < 6; i++)
    (*lstream)<<K_face[i] / face_area[i]<<" ";
  (*lstream)<<"\n";

  
  (*lstream)<<"12 edges lengthes are ";
  for(int i = 0; i < 12; i++)
    (*lstream)<<edge_length[i]<<" ";
  (*lstream)<<"\n";

  (*lstream)<<"K on 12 edges lengthes are ";
  for(int i = 0; i < 12; i++)
    (*lstream)<<K_edge[i] / edge_length[i]<<" ";
  (*lstream)<<"\n";

  
  (*lstream)<<"**********Finish outputing expansion info****** " <<"\n";
  
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
