#include "horizon.h"
#include "../../cosmo_includes.h"
#include "../../utils/math.h"
#include "../bssn/bssn_data.h"
#include <complex.h>
#include "../../utils/Eigen/Dense"

using namespace SAMRAI;

namespace cosmo
{
  
Horizon::Horizon(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  const tbox::Dimension& dim_in,
  boost::shared_ptr<tbox::Database> database_in,
  std::ostream* l_stream_in,
  int w_idx_in):
  lstream(l_stream_in),
  cosmo_horizon_db(database_in),
  dim(dim_in),
  is_periodic(cosmo_horizon_db->getBoolWithDefault("is_periodic", false)),
  is_sphere(cosmo_horizon_db->getBoolWithDefault("is_sphere", false)),
  radius_limit(cosmo_horizon_db->getDoubleWithDefault("radius_limit", INF)),
  w_idx(w_idx_in),
  invalid_id( hier::LocalId::getInvalidId(), tbox::SAMRAI_MPI::getInvalidRank())
{
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

  boost::shared_ptr<hier::VariableContext> context_scratch(
    variable_db->getContext("SCRATCH"));
  boost::shared_ptr<hier::VariableContext> context_active(
    variable_db->getContext("ACTIVE"));
  boost::shared_ptr<hier::VariableContext> context_previous(
    variable_db->getContext("PREVIOUS"));
  boost::shared_ptr<hier::VariableContext> context_k1(
    variable_db->getContext("RK_K1"));
  boost::shared_ptr<hier::VariableContext> context_k2(
    variable_db->getContext("RK_K2"));
  boost::shared_ptr<hier::VariableContext> context_k3(
    variable_db->getContext("RK_K3"));
  boost::shared_ptr<hier::VariableContext> context_k4(
    variable_db->getContext("RK_K4"));

  
  VAR_INIT(F);
  VAR_INIT(s1);
  VAR_INIT(s2);
  VAR_INIT(s3);

  
  REG_TO_CONTEXT(F, context_scratch, s, STENCIL_ORDER);
  REG_TO_CONTEXT(F, context_previous, p, STENCIL_ORDER);
  REG_TO_CONTEXT(F, context_active, a, STENCIL_ORDER);
  REG_TO_CONTEXT(F, context_k1, k1, STENCIL_ORDER);
  REG_TO_CONTEXT(F, context_k2, k2, STENCIL_ORDER);
  REG_TO_CONTEXT(F, context_k3, k3, STENCIL_ORDER);
  REG_TO_CONTEXT(F, context_k4, k4, STENCIL_ORDER);

  REG_TO_CONTEXT(s1, context_active, a, STENCIL_ORDER);
  REG_TO_CONTEXT(s2, context_active, a, STENCIL_ORDER);
  REG_TO_CONTEXT(s3, context_active, a, STENCIL_ORDER);

  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;


  const double * lower = &grid_geometry.getXLower()[0];
  const double * upper = &grid_geometry.getXUpper()[0];

  
  for(int i = 0 ; i < DIM; i++)
  {
    domain_lower[i] = lower[i];
    domain_upper[i] = upper[i];
  }
  origin = cosmo_horizon_db->getDoubleVector("origin");
  radius = cosmo_horizon_db->getDouble("radius");

  
  patch_work_i = patch_work_j = patch_work_k = -1;
  
}
Horizon::~Horizon()
{
  
}

void Horizon::clear(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln)
{
  math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);
  RK4_ARRAY_ZERO(F, hcellmath);
  EXTRA_ARRAY_ZERO(s2, hcellmath);
  EXTRA_ARRAY_ZERO(s2, hcellmath);
  EXTRA_ARRAY_ZERO(s3, hcellmath);
}

bool Horizon::initSphericalSurface(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  BSSN * bssn, boost::shared_ptr<hier::RefineOperator> space_refine_op)
{
  if(!is_sphere)
    TBOX_ERROR("The shape of apparent horizon is not set to sphere!\n");

  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());

  const_radius = 0;
  int cnt = 0;
  
  addNormVector(bssn, hierarchy);

  bool has_surface = false;
  
  int ln_num = hierarchy->getNumberOfLevels();

  double global_min_radius = INF;
  
  for(int ln = 0 ; ln < ln_num; ln ++)
  {
    boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

    for (hier::PatchLevel::iterator p(level->begin());
         p != level->end(); ++p)
    {
      const boost::shared_ptr<hier::Patch>& patch = *p;
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));
  
      initPData(patch);
      initMDA(patch);
      bssn->initPData(patch);
      bssn->initMDA(patch);

      boost::shared_ptr<pdat::CellData<real_t>> w_pdata (    
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  
          patch->getPatchData(w_idx)));
      arr_t w = pdat::ArrayDataAccess::access<DIM, double>(  
        w_pdata->getArrayData());
    
      const hier::Box& box = patch->getBox();
  
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];

      const real_t * dx = &(patch_geometry->getDx())[0];

#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSNData bd = {0};

            bssn->set_bd_values(i, j, k, &bd, dx);
            HorizonData hd = getHorizonData(i, j, k, &bd, dx);

            F_p(i, j, k) =
              ev_F(&bd, &hd, dx) / (bd.chi *
    sqrt(hd.d1F * hd.d1F * bd.gammai11 + hd.d2F * hd.d2F * bd.gammai22 + hd.d3F * hd.d3F * bd.gammai33
         + 2.0 * (hd.d1F * hd.d2F * bd.gammai12 + hd.d1F * hd.d3F * bd.gammai13 + hd.d2F * hd.d3F * bd.gammai23 )));

          }
        }
      }

#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            F_a(i, j, k) = F_p(i, j, k);
          }
        }
      }

    }

    xfer::RefineAlgorithm refiner;
    boost::shared_ptr<xfer::RefineSchedule> refine_schedule;

  
    registerRKRefinerActive(refiner, space_refine_op);

    level->getBoxLevel()->getMPI().Barrier();
    refine_schedule =
      refiner.createSchedule(level,
                             level,
                             ln-1,
                             hierarchy,
                             NULL);

    
    level->getBoxLevel()->getMPI().Barrier();
    refine_schedule->fillData(0.0);
    level->getBoxLevel()->getMPI().Barrier();

    for (hier::PatchLevel::iterator p(level->begin());
         p != level->end(); ++p)
    {
      const boost::shared_ptr<hier::Patch>& patch = *p;
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));
  
      initPData(patch);
      initMDA(patch);
      bssn->initPData(patch);
      bssn->initMDA(patch);

      boost::shared_ptr<pdat::CellData<real_t>> w_pdata (    
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  
          patch->getPatchData(w_idx)));
      arr_t w = pdat::ArrayDataAccess::access<DIM, double>(  
        w_pdata->getArrayData());
    
      const hier::Box& box = patch->getBox();
  
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];

      const real_t * dx = &(patch_geometry->getDx())[0];

#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            if(w(i, j, k) > 0 && onTheSurface(i, j, k))
            {
              real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0;
              real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0;
              real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0;
              real_t min_r = sqrt(pw2(x - origin[0]) + pw2(y - origin[1]) + pw2(z - origin[2]));

              if(min_r < radius_limit && min_r < global_min_radius)
                global_min_radius = min_r;
            }
          }
        }
      }
    }

       
  }

  mpi.Barrier();
  if (mpi.getSize() > 1) {
    mpi.AllReduce(&global_min_radius, 1, MPI_MIN);
  }
  mpi.Barrier();

  
  for(int ln = 0 ; ln < ln_num; ln ++)
  {
    boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

    for (hier::PatchLevel::iterator p(level->begin());
         p != level->end(); ++p)
    {
      const boost::shared_ptr<hier::Patch>& patch = *p;
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));
  
      initPData(patch);
      initMDA(patch);
      bssn->initPData(patch);
      bssn->initMDA(patch);

      boost::shared_ptr<pdat::CellData<real_t>> w_pdata (    
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  
          patch->getPatchData(w_idx)));
      arr_t w = pdat::ArrayDataAccess::access<DIM, double>(  
        w_pdata->getArrayData());
    
      const hier::Box& box = patch->getBox();
  
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];

      const real_t * dx = &(patch_geometry->getDx())[0];

      double sum_radius = 0;

      
#pragma omp parallel for collapse(2) reduction(+:sum_radius, cnt)      
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            if(w(i, j, k) > 0 && onTheSurface(i, j, k))
            {
              real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0;
              real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0;
              real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0;
              real_t min_r = sqrt(pw2(x - origin[0]) + pw2(y - origin[1]) + pw2(z - origin[2]));
              if(min_r < global_min_radius * 1.2 && min_r > global_min_radius * 0.8)
              {
                sum_radius += min_r;
                cnt ++;
              }
            }

          }
        }
      }

      const_radius += sum_radius;
    }    

  }

  mpi.Barrier();
  if (mpi.getSize() > 1) {
    mpi.AllReduce(&const_radius, 1, MPI_SUM);
    mpi.AllReduce(&cnt, 1, MPI_SUM);
  }
  mpi.Barrier();
  const_radius = const_radius / (double) cnt;

  if(cnt > 0) has_surface = true;
  return has_surface;
  
}
  
void Horizon::initSurface(const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{

  int ln_num = hierarchy->getNumberOfLevels();
  if(is_periodic)
  {
    for (int ln = 0; ln < ln_num; ln++)
    {
      boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    
      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
        const boost::shared_ptr<hier::Patch>& patch = *p;
        boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
          BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch->getPatchGeometry()));


        initPData(patch);
        initMDA(patch);
        const hier::Box& box = F_a_pdata->getGhostBox();
  
        const int * lower = &box.lower()[0];
        const int * upper = &box.upper()[0];

        const real_t * dx = &(patch_geometry->getDx())[0];

     
#pragma omp parallel for collapse(2)        
        for(int k = lower[2]; k <= upper[2]; k++)
        {
          for(int j = lower[1]; j <= upper[1]; j++)
          {
            for(int i = lower[0]; i <= upper[0]; i++)
            {
              real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0;
              real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0;
              real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0;
              real_t min_r = pw2(x - origin[0]) + pw2(y - origin[1]) + pw2(z - origin[2]);
              // min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_lower[0]) + pw2(y - domain_lower[1]) + pw2(z - domain_lower[2]));
              // min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_upper[0]) + pw2(y - domain_lower[1]) + pw2(z - domain_lower[2]));
              // min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_lower[0]) + pw2(y - domain_upper[1]) + pw2(z - domain_lower[2]));
              // min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_lower[0]) + pw2(y - domain_lower[1]) + pw2(z - domain_upper[2]));
              // min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_lower[0]) + pw2(y - domain_upper[1]) + pw2(z - domain_upper[2]));
              // min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_upper[0]) + pw2(y - domain_lower[1]) + pw2(z - domain_upper[2]));
              // min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_upper[0]) + pw2(y - domain_upper[1]) + pw2(z - domain_lower[2]));
              // min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_upper[0]) + pw2(y - domain_upper[1]) + pw2(z - domain_upper[2]));
              min_r = sqrt(min_r);
            
              F_a(i, j, k) = F_p(i, j, k) =
                min_r - radius; 
            }
          }
        }     
      }

    }  // loop over levels

    return;
  }
  for (int ln = 0; ln < ln_num; ln++)
  {
    boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    
    for (hier::PatchLevel::iterator p(level->begin());
         p != level->end(); ++p) {
      const boost::shared_ptr<hier::Patch>& patch = *p;
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));


      initPData(patch);
      initMDA(patch);
      const hier::Box& box = F_a_pdata->getGhostBox();
  
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];

      const real_t * dx = &(patch_geometry->getDx())[0];

     
      #pragma omp parallel for collapse(2)              
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0;
            real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0;
            real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0;
            F_a(i, j, k) = F_p(i, j, k) =
              sqrt(pw2(x - origin[0]) + pw2(y - origin[1]) + pw2(z - origin[2])) - radius;
          }
        }
      }     
    }

  }  // loop over levels

}

/**
 * @brief initilizing all pointers for every components's patchdata
 * 
 */
void Horizon::initPData(
  const boost::shared_ptr<hier::Patch> & patch)
{
  PDATA_ALL_INIT(F);
  PDATA_INIT(s1, a);
  PDATA_INIT(s2, a);
  PDATA_INIT(s3, a);
}


/**
 * @brief initilizing all pointers for every components's array access
 * 
 */
void Horizon::initMDA(
  const boost::shared_ptr<hier::Patch> & patch)
{

  MDA_ACCESS_ALL_INIT(F);
  MDA_ACCESS_INIT(s1, a);
  MDA_ACCESS_INIT(s2, a);
  MDA_ACCESS_INIT(s3, a);
}

void Horizon::setLevelTime(
  const boost::shared_ptr<hier::PatchLevel> & level,
  double from_t, double to_t)
{
  SET_LEVEL_TIME(F, from_t, to_t);
}

  
HorizonData Horizon::getHorizonData(
  idx_t i, idx_t j, idx_t k, BSSNData *bd, const real_t dx[])
{
  HorizonData hd = {0};

  hd.F = F_a(i, j, k);
  hd.s1 = s1_a(i, j, k);
  hd.s2 = s2_a(i, j, k);
  hd.s3 = s3_a(i, j, k);

  hd.d1F = derivative(i, j, k, 1, F_a, dx);
  hd.d2F = derivative(i, j, k, 2, F_a, dx);
  hd.d3F = derivative(i, j, k, 3, F_a, dx);

  hd.d1d1F = double_derivative(i, j, k, 1, 1, F_a, dx);
  hd.d2d2F = double_derivative(i, j, k, 2, 2, F_a, dx);
  hd.d3d3F = double_derivative(i, j, k, 3, 3, F_a, dx);
  hd.d1d2F = double_derivative(i, j, k, 1, 2, F_a, dx);
  hd.d1d3F = double_derivative(i, j, k, 1, 3, F_a, dx);
  hd.d2d3F = double_derivative(i, j, k, 2, 3, F_a, dx);
  
  return hd;
}

void Horizon::addNormVector(
  BSSN * bssn, const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  for(idx_t ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
  {
    const boost::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));
    addNormVector(bssn, level);
  }
  
}

void Horizon::addNormVector(
  BSSN * bssn, const boost::shared_ptr<hier::PatchLevel> & level)
{
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    addNormVector(bssn, patch);
  }  
}

  
void Horizon::addNormVector(
  BSSN * bssn, const boost::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);
  bssn->initPData(patch);
  bssn->initMDA(patch);

  
  const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom( 
    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  const real_t * dx = &(patch_geom->getDx())[0];

  const hier::Box& box = patch->getBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

      #pragma omp parallel for collapse(2)          
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        BSSNData bd = {0};
        // TODO: remove redundant computations here?
        bssn->set_bd_values(i, j, k, &bd, dx);
        HorizonData hd = getHorizonData(i, j, k, &bd, dx);

        s1_a(i, j, k) = bd.chi *
          (bd.gammai11 * hd.d1F + bd.gammai12 * hd.d2F + bd.gammai13 * hd.d3F) /
          sqrt( 
               (bd.gammai11 * hd.d1F * hd.d1F + bd.gammai22 * hd.d2F * hd.d2F + bd.gammai33 * hd.d3F *hd.d3F
                + 2.0 * (bd.gammai12 * hd.d1F * hd.d2F + bd.gammai13 * hd.d1F * hd.d3F + bd.gammai23 * hd.d2F * hd.d3F)));

        s2_a(i, j, k) = bd.chi *
          (bd.gammai21 * hd.d1F + bd.gammai22 * hd.d2F + bd.gammai23 * hd.d3F) /
          sqrt(
               (bd.gammai11 * hd.d1F * hd.d1F + bd.gammai22 * hd.d2F * hd.d2F + bd.gammai33 * hd.d3F *hd.d3F
                + 2.0 * (bd.gammai12 * hd.d1F * hd.d2F + bd.gammai13 * hd.d1F * hd.d3F + bd.gammai23 * hd.d2F * hd.d3F)));
        s3_a(i, j, k) = bd.chi *
          (bd.gammai31 * hd.d1F + bd.gammai32 * hd.d2F + bd.gammai33 * hd.d3F) /
          sqrt(
               (bd.gammai11 * hd.d1F * hd.d1F + bd.gammai22 * hd.d2F * hd.d2F + bd.gammai33 * hd.d3F *hd.d3F
                + 2.0 * (bd.gammai12 * hd.d1F * hd.d2F + bd.gammai13 * hd.d1F * hd.d3F + bd.gammai23 * hd.d2F * hd.d3F)));

      }
    }
  }
  return;      
}

real_t Horizon::ev_F(BSSNData *bd, HorizonData *hd, const real_t dx[])
{
  real_t diFdiF = pw2(bd->chi) * (hd->d1F * hd->d1F * bd->gammai11 + hd->d2F * hd->d2F * bd->gammai22 + hd->d3F * hd->d3F * bd->gammai33
         + 2.0 * (hd->d1F * hd->d2F * bd->gammai12 + hd->d1F * hd->d3F * bd->gammai13 + hd->d2F * hd->d3F * bd->gammai23 ));

  real_t kappa =  ((pw2(bd->chi) * COSMO_SUMMATION_2(HORIZON_CALCULATE_DS0)) / sqrt(diFdiF)
           - pw2(bd->chi) * pw2(bd->chi) * COSMO_SUMMATION_2(HORIZON_CALCULATE_DS) / pow(diFdiF, 1.5)
    - bd->K
    +  1.0 / pw2(bd->chi) * (
      (bd->A11 + bd->K * bd->gamma11 / 3.0) * hd->s1 * hd->s1
      + (bd->A22 + bd->K * bd->gamma22 / 3.0) * hd->s2 * hd->s2
      + (bd->A33 + bd->K * bd->gamma33 / 3.0) * hd->s3 * hd->s3
      + 2.0 * (bd->A12 + bd->K * bd->gamma12 / 3.0) * hd->s1 * hd->s2
      + 2.0 * (bd->A13 + bd->K * bd->gamma13 / 3.0) * hd->s1 * hd->s3
      + 2.0 * (bd->A23 + bd->K * bd->gamma23 / 3.0) * hd->s2 * hd->s3
    )) * bd->chi *
    sqrt(hd->d1F * hd->d1F * bd->gammai11 + hd->d2F * hd->d2F * bd->gammai22 + hd->d3F * hd->d3F * bd->gammai33
         + 2.0 * (hd->d1F * hd->d2F * bd->gammai12 + hd->d1F * hd->d3F * bd->gammai13 + hd->d2F * hd->d3F * bd->gammai23 ));

  real_t x = domain_lower[0] + (double)bd->i * dx[0] + dx[0]/2.0;
  real_t y = domain_lower[1] + (double)bd->j * dx[1] + dx[1]/2.0;
  real_t z = domain_lower[2] + (double)bd->k * dx[2] + dx[2]/2.0;

  if(bd->i == 105 && bd->j == 64 && bd->k ==64)
  {
    std::cout<<kappa <<"\n";
  }
                
  
  return kappa ;
}
void Horizon::RKEvolvePt(
  idx_t i, idx_t j, idx_t k, BSSNData &bd, const real_t dx[], real_t dt)
{
  HorizonData hd = getHorizonData(i, j, k, &bd, dx);

  F_s(i,j,k) = ev_F(&bd, &hd, dx) * dt;
}

bool Horizon::onTheSurface(idx_t i, idx_t j, idx_t k)
{
  if(SIGN(F_a(i-1, j, k)) * SIGN(F_a(i+1, j, k)) < 0 
     || SIGN(F_a(i, j-1, k)) * SIGN(F_a(i, j+1, k)) < 0
     || SIGN(F_a(i, j, k-1)) * SIGN(F_a(i, j, k+1)) < 0)
    return true;
  
  return false;
}

bool Horizon::belowTheSurface(idx_t i, idx_t j, idx_t k)
{
  if(F_a(i, j, k) < 0 &&
     (F_a(i+1, j, k) > 0 || F_a(i-1, j, k) > 0 ||
      F_a(i, j+1, k) > 0 || F_a(i, j-1, k) > 0 ||
      F_a(i, j, k+1) > 0 || F_a(i, j, k-1) > 0))
    return true;
  
  return false;
}

void Horizon::prepareForK1(
  const boost::shared_ptr<hier::PatchLevel> & level,
  real_t to_t)
{
  if(level == NULL) return;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    
    initPData(patch);
    initMDA(patch);

    const hier::Box& box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];


    // right branch of the tree
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(F_a_idx)->getTime()) < EPS)
    {
            #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_R_K1(F);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(F_p_idx)->getTime())
       - (patch->getPatchData(F_a_idx)->getTime()
          - patch->getPatchData(F_p_idx)->getTime())) < EPS)
    {
            #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_L_K1(F);
          }
        }
      }
    }
    else
      TBOX_ERROR("Current level locates neigher L or R branch of its father, check your code!");
    
  }
}

/**
 * @brief calculate K2 value on coarser level
 * 
 */
void Horizon::prepareForK2(
  const boost::shared_ptr<hier::PatchLevel> & level,
  real_t to_t)
{
  if(level == NULL) return;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    initPData(patch);
    initMDA(patch);

    const hier::Box& box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];


    // right branch of the tree
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(F_a_idx)->getTime()) < EPS)
    {
#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_R_K2(F);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(F_p_idx)->getTime())
       - (patch->getPatchData(F_a_idx)->getTime()
          - patch->getPatchData(F_p_idx)->getTime())) < EPS)
    {
#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_L_K2(F);
          }
        }
      }
    }
  }
}

/**
 * @brief calculate K3 value on coarser level
 * 
 */  
void Horizon::prepareForK3(
  const boost::shared_ptr<hier::PatchLevel> & level,
  real_t to_t)
{
  if(level == NULL) return;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    initPData(patch);
    initMDA(patch);

    const hier::Box& box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    // right branch of the tree
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(F_a_idx)->getTime()) < EPS)
    {
#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_R_K3(F);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(F_p_idx)->getTime())
       - (patch->getPatchData(F_a_idx)->getTime()
          - patch->getPatchData(F_p_idx)->getTime())) < EPS)
    {
#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_L_K3(F);
          }
        }
      }
    }
  }
}

/**
 * @brief calculate K4 value on coarser level
 * 
 */
void Horizon::prepareForK4(
  const boost::shared_ptr<hier::PatchLevel> & level,
  double to_t)
{
  if(level == NULL) return;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    initPData(patch);
    initMDA(patch);

    const hier::Box& box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    
    


    // right branch of the tree
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(F_a_idx)->getTime()) < EPS)
    {
            #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_R_K4(F);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(F_p_idx)->getTime())
       - (patch->getPatchData(F_a_idx)->getTime()
          - patch->getPatchData(F_p_idx)->getTime())) < EPS)
    {
#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_L_K4(F);
          }
        }
      }
    }
  }
}  

void Horizon::alloc(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln)
{
  boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
  RK4_ARRAY_ALLOC(F);
  EXTRA_ARRAY_ALLOC(s1);
  EXTRA_ARRAY_ALLOC(s2);
  EXTRA_ARRAY_ALLOC(s3);
}


void Horizon::registerRKRefinerActive(
  xfer::RefineAlgorithm& refiner,
  boost::shared_ptr<hier::RefineOperator> &space_refine_op)
{
  REGISTER_SPACE_REFINE_A(F, refiner, space_refine_op);
}

  
/**
 * @brief register "scratch" component of the fields to refiner
 * 
 * @param refiner
 * @param refiner operator
 */
void Horizon::registerRKRefiner(
  xfer::RefineAlgorithm& refiner,
  boost::shared_ptr<hier::RefineOperator> &space_refine_op)
{
  REGISTER_SPACE_REFINE_S(F, refiner, space_refine_op);
}

/**
 * @brief register "active" component of the fields to coarsener 
 * 
 * @param coarsener
 * @param coarse operator
 */
void Horizon::registerCoarsenActive(
  xfer::CoarsenAlgorithm& coarsener,
  boost::shared_ptr<hier::CoarsenOperator>& coarsen_op)
{
  REGISTER_COARSEN_A(F, coarsener, coarsen_op);
}


void Horizon::copyAToP(
  math::HierarchyCellDataOpsReal<real_t> & hcellmath)
{
  COPY_A_TO_P(F);
}

/**
 * @brief  finalize k1 step for RK on both interior and boundary
 * 
 */
void Horizon::K1FinalizePatch(
  const boost::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);
  
  const hier::Box& box = F_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

      #pragma omp parallel for collapse(2)          
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        RK4_FINALIZE_FIELD_1(F);
      }
    }
  }
}

void Horizon::K2FinalizePatch(
  const boost::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  const hier::Box& box = F_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

      #pragma omp parallel for collapse(2)          
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        RK4_FINALIZE_FIELD_2(F);
      }
    }
  }
}


void Horizon::K3FinalizePatch(
  const boost::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  const hier::Box& box = F_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

      #pragma omp parallel for collapse(2)          
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        RK4_FINALIZE_FIELD_3(F);
      }
    }
  }
}


void Horizon::K4FinalizePatch(
  const boost::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  //Only evolv inner grids for K4
  const hier::Box& box = F_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

        #pragma omp parallel for collapse(2)        
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        RK4_FINALIZE_FIELD_4(F);
      }
    }
  }
}


void Horizon::RKEvolveHorizon(
  const boost::shared_ptr<hier::Patch> & patch, BSSN * bssn, real_t dt)
{
  bssn->initPData(patch);
  bssn->initMDA(patch);
  initPData(patch);
  initMDA(patch);

  
  const hier::Box& box = patch->getBox();
  
  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
  const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(  
    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  //initialize dx for each patch
  const real_t * dx = &(patch_geom->getDx())[0];


#pragma omp parallel for collapse(2)          
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        if( i == 105 && j == 64 && k == 64)
          std::cout<<"Level set function at closest point is "<<F_a(i, j, k)<<"\n";
        BSSNData bd = {0};

        bssn->set_bd_values(i, j, k, &bd, dx);
        RKEvolvePt(i, j, k, bd, dx, dt);
        
      }
    }
  }
}

real_t Horizon::maxSurfaceMove(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln, idx_t w_idx)
{
  real_t max_change = -1;
  boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

  
  for (hier::PatchLevel::iterator p(level->begin());
       p != level->end(); ++p)
  {
    const boost::shared_ptr<hier::Patch>& patch = *p;
    boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
      BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry()));
  
    initPData(patch);

    initMDA(patch);

    boost::shared_ptr<pdat::CellData<real_t>> w_pdata (    
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  
        patch->getPatchData(w_idx)));
    arr_t w = pdat::ArrayDataAccess::access<DIM, double>(  
      w_pdata->getArrayData());
    
    const hier::Box& box = patch->getBox();
  
    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    const real_t * dx = &(patch_geometry->getDx())[0];
    
#pragma omp parallel for collapse(2) reduction(max : max_change)
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          if(w(i, j, k) > 0 && onTheSurface(i, j, k)
             && fabs(F_a(i, j, k) - F_p(i, j, k)) > max_change)
            max_change = fabs(F_a(i, j, k) - F_p(i, j, k));
        }
      }
    }


  }
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  mpi.Barrier();
  if (mpi.getSize() > 1) {
    mpi.AllReduce(&max_change, 1, MPI_MAX);
  }

  return max_change;
  
}

real_t Horizon::getFonBD(
  const boost::shared_ptr<hier::Patch> & patch,
  real_t bx, real_t by, real_t bz,
  hier::Box& g_box)
{
  initPData(patch);
  initMDA(patch);

  const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom( 
    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  const real_t * dx = &(patch_geom->getDx())[0];

  
  const int * lower = &g_box.lower()[0];
  const int * upper = &g_box.upper()[0];

  real_t min_s = INF;
  idx_t min_i = lower[0], min_j = lower[1], min_k = lower[2];
  real_t min_r = -1.0;
  
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0;
        real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0;
        real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0;

        real_t a = (x - origin[0]) * (by - origin[1]) - (y - origin[1]) * (bx - origin[0]);
        real_t b = (x - origin[0]) * (bz - origin[2]) - (z - origin[2]) * (bx - origin[0]);
        real_t c = (y - origin[1]) * (bz - origin[2]) - (z - origin[2]) * (by - origin[1]);
        
        real_t  slope = tbox::MathUtilities<double>::Max(
          tbox::MathUtilities<double>::Max(
            tbox::MathUtilities<double>::Abs(a), tbox::MathUtilities<double>::Abs(b)),
          tbox::MathUtilities<double>::Abs(c));
        if(slope < min_s)
        {
          min_s = slope;
          min_i = i, min_j = j, min_k = k;
          min_r = sqrt(pw2(x - origin[0]) + pw2(y - origin[1]) + pw2(z - origin[2]));
        }
      }
    }
  }

  if(min_r < 0)
    TBOX_ERROR("Cannot find grid with right slope!");

  return min_r - F_a(min_i, min_j, min_k);
  
}
  
void Horizon::updateBD(
  const boost::shared_ptr<hier::Patch> & patch,
  real_t dt)
{
  boost::shared_ptr<hier::PatchGeometry> geom (patch->getPatchGeometry());

  idx_t codim = 1;

  // getting all codimension 1 boxes
  const std::vector<hier::BoundaryBox> & codim1_boxes =
    geom->getCodimensionBoundaries(codim);

  const idx_t n_codim1_boxes = static_cast<idx_t>(codim1_boxes.size());

  // if it has no codimention 1 boundary, it has no other type of boundaries
  if(n_codim1_boxes == 0) return;

  initPData(patch);
  initMDA(patch);

  hier::Box& gg_box = const_cast<hier::Box&> (F_a_pdata->getGhostBox());

  hier::Box g_box = gg_box;
  
  const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom( 
    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

    
  const hier::Box& patch_box = patch->getBox();

  for(int i = 0 ; i < n_codim1_boxes; i++)
  {
    hier::Box boundary_fill_box =
      geom->getBoundaryFillBox(
        codim1_boxes[i], patch_box, F_a_pdata->getGhostCellWidth());

    if(boundary_fill_box.empty()) continue;


    idx_t l_idx = codim1_boxes[i].getLocationIndex();

    if(l_idx%2)
      g_box.growUpper(
        (hier::Box::dir_t)l_idx/2, -2*STENCIL_ORDER_WIDTH);
    else
      g_box.growLower(
        (hier::Box::dir_t)l_idx/2, -2*STENCIL_ORDER_WIDTH);
  }

  
  for(int l = 0 ; l < n_codim1_boxes; l++)
  {
    hier::Box boundary_fill_box =
      geom->getBoundaryFillBox(
        codim1_boxes[l], patch_box, F_a_pdata->getGhostCellWidth());

    if(boundary_fill_box.empty()) continue;


    idx_t l_idx = codim1_boxes[l].getLocationIndex();
    
    BSSNData bd = {0};

    boundary_fill_box.shift(
      (hier::Box::dir_t)l_idx/2,
      (l_idx%2)?(-STENCIL_ORDER_WIDTH):STENCIL_ORDER_WIDTH);


    boundary_fill_box *= patch_box;

    //initialize dx for each patch
    const real_t * dx = &(patch_geom->getDx())[0];

    const idx_t * lower = &boundary_fill_box.lower()[0];
    const idx_t * upper = &boundary_fill_box.upper()[0];
    if(l_idx == 0)
    {
#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = upper[0]; i >= lower[0]; i--)
          {
            F_a(i, j, k) = 3.0 * (F_a(i+1, j, k) - F_a(i+2, j, k)) + F_a(i+3, j, k);

          }
        }
      }
    }
    if(l_idx == 1)
    {
#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            F_a(i, j, k) = 3.0 * (F_a(i-1, j, k) - F_a(i-2, j, k)) + F_a(i-3, j, k);

          }
        }
      }
    }
    if(l_idx == 2)
    {
#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = upper[1]; j >= lower[1]; j--)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            F_a(i, j, k) = 3.0 * (F_a(i, j+1, k) - F_a(i, j+2, k)) + F_a(i, j+3, k);

          }
        }
      }
    }
    if(l_idx == 3)
    {
#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            F_a(i, j, k) = 3.0 * (F_a(i, j-1, k) - F_a(i, j-2, k)) + F_a(i, j-3, k);

          }
        }
      }
    }
    if(l_idx == 4)
    {
#pragma omp parallel for collapse(2)        
      for(int k = upper[2]; k >= lower[2]; k--)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            F_a(i, j, k) = 3.0 * (F_a(i, j, k+1) - F_a(i, j, k+2)) + F_a(i, j, k+3);

          }
        }
      }
    }
    if(l_idx == 5)
    {
#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            F_a(i, j, k) = 3.0 * (F_a(i, j, k-1) - F_a(i, j, k-2)) + F_a(i, j, k-3);

          }
        }
      }
    }

    
  }

//   /************************updating codim = 2 boundaries****************/
//   codim = 2;

//   const std::vector<hier::BoundaryBox> & codim2_boxes =
//     geom->getCodimensionBoundaries(codim);

//   const idx_t n_codim2_boxes = static_cast<idx_t>(codim2_boxes.size());




//   for(int l = 0 ; l < n_codim2_boxes; l++)
//   {
//     hier::Box  boundary_fill_box =
//       geom->getBoundaryFillBox(
//         codim2_boxes[l], patch_box, F_a_pdata->getGhostCellWidth());
    
//     if(boundary_fill_box.empty()) continue;  


//     idx_t l_idx = codim2_boxes[l].getLocationIndex();
    
//     BSSNData bd = {0};

//     std::vector<idx_t> shift_vec;
//     if(l_idx == 0 || l_idx == 2 || l_idx == 4 || l_idx == 6)
//       shift_vec.push_back(STENCIL_ORDER_WIDTH);
//     else if(l_idx == 1 || l_idx == 3 || l_idx == 5 || l_idx == 7)
//       shift_vec.push_back(-STENCIL_ORDER_WIDTH);
//     else
//       shift_vec.push_back(0);
    
//     if(l_idx == 0 || l_idx == 1 || l_idx == 8 || l_idx == 10)
//       shift_vec.push_back(STENCIL_ORDER_WIDTH);
//     else if(l_idx == 2 || l_idx == 3 || l_idx ==9 || l_idx == 11)
//      shift_vec.push_back(-STENCIL_ORDER_WIDTH);
//     else
//       shift_vec.push_back(0);

//     if( l_idx == 4 || l_idx == 5 || l_idx == 8 || l_idx == 9)
//       shift_vec.push_back(STENCIL_ORDER_WIDTH);
//     else if(l_idx == 6 || l_idx == 7 || l_idx == 10 || l_idx == 11)
//       shift_vec.push_back(-STENCIL_ORDER_WIDTH);
//     else
//       shift_vec.push_back(0);

//     boundary_fill_box.shift(hier::IntVector(shift_vec));

//     boundary_fill_box *= patch_box;

//     //initialize dx for each patch
//     const real_t * dx = &(patch_geom->getDx())[0];

//     const idx_t * lower = &boundary_fill_box.lower()[0];
//     const idx_t * upper = &boundary_fill_box.upper()[0];

    
//     for(int k = lower[2]; k <= upper[2]; k++)
//     {
//       for(int j = lower[1]; j <= upper[1]; j++)
//       {
//         for(int i = lower[0]; i <= upper[0]; i++)
//         {
//           real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0;
//           real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0;
//           real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0;
        
//           F_a(i, j, k) = sqrt(pw2(x - origin[0]) + pw2(y - origin[1]) + pw2(z - origin[2]))
//             - getFonBD(patch, x, y, z, g_box);
//         }
//       }
//     }
//   }

//   /************************updating codim = 3 boundaries****************/
//   codim = 3;

//   const std::vector<hier::BoundaryBox> & codim3_boxes =
//     geom->getCodimensionBoundaries(codim);

//   const idx_t n_codim3_boxes = static_cast<idx_t>(codim3_boxes.size());




//   for(int l = 0 ; l < n_codim3_boxes; l++)
//   {
//     hier::Box boundary_fill_box =
//       geom->getBoundaryFillBox(
//         codim3_boxes[l], patch_box, F_a_pdata->getGhostCellWidth());

//     if(boundary_fill_box.empty()) continue;

//     idx_t l_idx = codim3_boxes[l].getLocationIndex();

//     std::vector<idx_t> shift_vec;
    
//     if(l_idx == 0 || l_idx == 2 || l_idx == 4 || l_idx == 6)
//       shift_vec.push_back(STENCIL_ORDER_WIDTH);
//     else
//       shift_vec.push_back(-STENCIL_ORDER_WIDTH);
    
//     if(l_idx == 0 || l_idx == 1 || l_idx == 4 || l_idx == 5)
//       shift_vec.push_back(STENCIL_ORDER_WIDTH);
//     else
//       shift_vec.push_back(-STENCIL_ORDER_WIDTH);

//     if( l_idx == 0 || l_idx == 1 || l_idx == 2 || l_idx == 3)
//       shift_vec.push_back(STENCIL_ORDER_WIDTH);
//     else
//       shift_vec.push_back(-STENCIL_ORDER_WIDTH);

//     boundary_fill_box.shift(hier::IntVector(shift_vec));

//     boundary_fill_box *= patch_box;

    
//     const idx_t * lower = &boundary_fill_box.lower()[0];
//     const idx_t * upper = &boundary_fill_box.upper()[0];

    
//     BSSNData bd = {0};


//     //initialize dx for each patch
//     const real_t * dx = &(patch_geom->getDx())[0];
  
//     for(int k = lower[2]; k <= upper[2]; k++)
//     {
//       for(int j = lower[1]; j <= upper[1]; j++)
//       {
//         for(int i = lower[0]; i <= upper[0]; i++)
//         {
//           real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0;
//           real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0;
//           real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0;
        
//           F_a(i, j, k) = sqrt(pw2(x - origin[0]) + pw2(y - origin[1]) + pw2(z - origin[2]))
//             - getFonBD(patch, x, y, z, g_box);
//         }
//       }
//     }
//   }

// }
}

// finding maximum radius of the horizon to set spacing of the mesh
// and initializing the maps from angle to radius
real_t Horizon::findMaxHorizonRadius(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double theta_0, double phi_0)
{
  if(is_sphere) return const_radius;
  real_t max_r = -1;
  int ln_num = hierarchy->getNumberOfLevels();

  for (int ln = 0; ln < ln_num; ln++)
  {
    boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    
    for (hier::PatchLevel::iterator p(level->begin());
         p != level->end(); ++p) {
      const boost::shared_ptr<hier::Patch>& patch = *p;
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));


      initPData(patch);
      initMDA(patch);
      const hier::Box& box = patch->getBox();
  
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];

      const real_t * dx = &(patch_geometry->getDx())[0];

      boost::shared_ptr<pdat::CellData<real_t>> w_pdata (    
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  
          patch->getPatchData(w_idx)));
      arr_t w = pdat::ArrayDataAccess::access<DIM, double>(  
        w_pdata->getArrayData());

      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            if(w(i, j, k) > 0 && belowTheSurface(i, j, k))
            {

              real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0 - origin[0];
              real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0 - origin[1];
              real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0 - origin[2];

              // distance plus correctino
              real_t r = sqrt(pw2(x ) + pw2(y ) + pw2(z )) - F_a(i, j, k);

              //              std::cout<<i<<" "<<j<<" "<<k<<" "<<r<<"\n";
              if(r < radius_limit)
              {
              
                max_r = tbox::MathUtilities<double>::Max(max_r, r);


              
                real_t phi = atan(y / x);
                real_t theta = acos(z / sqrt(pw2(x) + pw2(y) + pw2(z))); 

                //getting the minimum distance between grid cell and certain direction
                min_d = tbox::MathUtilities<double>::Min(
                  //min_d, sqrt(PW2(r *(theta - theta_0 ) ) + PW2(r*(phi-phi_0))));
                  min_d, r * acos(cos(theta) * cos(theta_0) + sin(theta) * sin(theta_0) * cos(phi-phi_0)));
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
    mpi.AllReduce(&max_r, 1, MPI_MAX);
    mpi.AllReduce(&min_d, 1, MPI_MIN);
  }
  mpi.Barrier();
  return max_r;
}



// try to find corresponding surface radius
// through interpolating from ordered radius map
// return -1 if no proper radius found (possible for some MPI threads)
real_t Horizon::findRadius(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double theta_0, double phi_0)
{

  if(is_sphere) return const_radius;
  real_t res = -1;
  min_d = INF;
  int ln_num = hierarchy->getNumberOfLevels();
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());

  // marks whether patch needs to be changed
  int need_to_change_patch = 0;

  if(mpi.getRank() == patch_work_mpi_rank)
  {
  for (int ln = 0; ln < ln_num; ln++)
  {
    boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    
    for (hier::PatchLevel::iterator p(level->begin());
         p != level->end(); ++p)
    {
      const boost::shared_ptr<hier::Patch>& patch = *p;
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      //if(patch->getGlobalId() != hier::GlobalId(hier::LocalId(local_id), patch_work_mpi_rank)) continue;

      if(patch->getLocalId() != hier::LocalId(local_id))
        continue;

      initPData(patch);
      initMDA(patch);

      const hier::Box& box = patch->getBox();
  
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];

      const real_t * dx = &(patch_geometry->getDx())[0];

      boost::shared_ptr<pdat::CellData<real_t>> w_pdata (    
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  
          patch->getPatchData(w_idx)));
      arr_t w = pdat::ArrayDataAccess::access<DIM, double>(  
        w_pdata->getArrayData());

      res = sqrt(pw2(domain_lower[0] + (double)patch_work_i * dx[0] + dx[0]/2.0 - origin[0] )
                 + pw2(domain_lower[1] + (double)patch_work_j * dx[1] + dx[1]/2.0 - origin[1] )
                 + pw2(domain_lower[2] + (double)patch_work_k * dx[2] + dx[2]/2.0 - origin[2] ))
      - F_a(patch_work_i, patch_work_j, patch_work_k);

      int lk = patch_work_k - STENCIL_ORDER+1, uk = patch_work_k + STENCIL_ORDER-1;
      int lj = patch_work_j - STENCIL_ORDER+1, uj = patch_work_j + STENCIL_ORDER-1;
      int li = patch_work_i - STENCIL_ORDER+1, ui = patch_work_i + STENCIL_ORDER-1;
      
      for(int k = lk; k <= uk; k++)
        for(int j = lj; j <= uj; j++)
          for(int i = li; i <= ui; i++)
          {

            if(w(i, j, k) > 0 && belowTheSurface(i, j, k))
            {
              real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0 - origin[0];
              real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0 - origin[1];
              real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0 - origin[2];

              real_t r = sqrt(pw2(x ) + pw2(y ) + pw2(z )) - F_a(i, j, k);
              real_t phi = atan(y / x);
              real_t theta = acos(z / sqrt(pw2(x) + pw2(y) + pw2(z))); 

              if(r * acos(cos(theta) * cos(theta_0) + sin(theta) * sin(theta_0) * cos(phi-phi_0)) < min_d)
              //if(sqrt(PW2(r *(theta - theta_0 ) ) + PW2(r*(phi-phi_0))) < min_d)
              {
                //min_d = sqrt(PW2(r *(theta - theta_0 ) ) + PW2(r*(phi-phi_0)));
                min_d = r * acos(cos(theta) * cos(theta_0) + sin(theta) * sin(theta_0) * cos(phi-phi_0));
                patch_work_i = i;
                patch_work_j = j;
                patch_work_k = k;
                res = r;
              }
            }
          }

   
      // it goes outside certain patch
      if(patch_work_i < lower[0] || patch_work_i > upper[0]
         || patch_work_j < lower[1] || patch_work_j > upper[1]
         || patch_work_k < lower[2] || patch_work_k > upper[2])
      {
        need_to_change_patch = 1;
      }


    }
  }
  }

  mpi.Barrier();

  if (mpi.getSize() > 1){

    mpi.Bcast(&min_d, 1, MPI_DOUBLE, patch_work_mpi_rank);

    mpi.Bcast(&patch_work_i, 1, MPI_INT, patch_work_mpi_rank);
    mpi.Bcast(&patch_work_j, 1, MPI_INT, patch_work_mpi_rank);
    mpi.Bcast(&patch_work_k, 1, MPI_INT, patch_work_mpi_rank);
    
    mpi.Bcast(&res, 1, MPI_DOUBLE, patch_work_mpi_rank);
    mpi.Bcast(&need_to_change_patch, 1, MPI_INT, patch_work_mpi_rank);

  }
  mpi.Barrier();

  
  if(need_to_change_patch)
  {

    findPatch(hierarchy, theta_0, phi_0);
    
  }

  return res;
  
}

  
void Horizon::set_G_values(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double theta, double phi, int theta_i, int phi_i, double r, KillingData *kd, BSSN * bssn)
{
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  cur_mpi_rank = -1;
  //if(patch_work_mpi_rank == mpi.getRank())
  int ln_num = hierarchy->getNumberOfLevels(), cur_mpi_level = -1, ln;

  real_t x = r * cos(phi) * sin(theta);
  real_t y = r * sin(phi) * sin(theta);
  real_t z = r * cos(theta);
        

  real_t st = sin(theta);
  real_t ct = cos(theta);
  real_t sp = sin(phi);
  real_t cp = cos(phi);

  for (ln = ln_num - 1; ln >= 0; ln--)
  {
    boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    hier::PatchLevel::iterator p(level->begin());
    for (;p != level->end(); ++p)
    {
      const boost::shared_ptr<hier::Patch>& patch = *p;

      const hier::Box& box = patch->getBox();

      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];
      
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));


      const real_t * dx = &(patch_geometry->getDx())[0];
      int i0 = floor((x + origin[0] - domain_lower[0] ) / dx[0] - 0.5);
      int j0 = floor((y + origin[1] - domain_lower[1] ) / dx[1] - 0.5);
      int k0 = floor((z + origin[2] - domain_lower[2] ) / dx[2] - 0.5);

      if( i0 >= lower[0] && i0 <= upper[0]
          && j0 >= lower[1] && j0 <= upper[1]
          && k0 >= lower[2] && k0 <= upper[2])
      {
        cur_mpi_rank = mpi.getRank();
        cur_mpi_level = ln;
        
        G111[theta_i][phi_i] = G112[theta_i][phi_i] = G122[theta_i][phi_i]
          = G211[theta_i][phi_i] = G212[theta_i][phi_i] = G222[theta_i][phi_i] = 0;


        
        
        initPData(patch);
        initMDA(patch);
        bssn->initPData(patch);
        bssn->initMDA(patch);
        
        BSSNData bd = {0};
        
        for(int k = k0; k <= k0 + 1; k++)
        {
          real_t z0 = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0 - origin[2];
          
          COSMO_APPLY_TO_IJK_PERMS(HORIZON_DEFINE_TEMP_GJ);
          
          for(int j = j0; j <= j0 + 1; j++)
          {
            real_t y0 = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0 - origin[1];
      
            COSMO_APPLY_TO_IJK_PERMS(HORIZON_DEFINE_TEMP_GI);

            for(int i = i0; i <= i0 + 1; i++)
            {
              real_t x0 = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0 - origin[0];
              bssn->set_bd_values(i, j, k, &bd, dx);
              COSMO_APPLY_TO_IJK_PERMS(HORIZON_INTERPOLATE_G_1);

            }
      
            COSMO_APPLY_TO_IJK_PERMS(HORIZON_INTERPOLATE_G_2);
          }
    
          COSMO_APPLY_TO_IJK_PERMS(HORIZON_INTERPOLATE_G_3);
        }
  
        G111[theta_i][phi_i] = ((pw2(x) + pw2(y))*z*(5*(pw2(x) + pw2(y)) - 
                                                     pw2(z))*cos(2*theta) + 
                                (pw2(x) + pw2(y))*(-(z*(3*(pw2(x) + pw2(y)) + 
                                                        pw2(z))) + 4*(cp*x + sp*y)*(pw2(x) + pw2(y) - 
                                                                                    pw2(z))*sin(2*theta)) + 
                                2*pw2(ct)*z*(3*(pw2(x) + pw2(y)) + pw2(z))*((x - 
                                                                             y)*(x + y)*cos(2*phi) + 2*x*y*sin(2*phi)) + 
                                pw2(r)*(pw2(x) + 
                                        pw2(y))*(4*pw2(cp)*pw2(ct)*x*z*kd->Gc111 + x*z*kd->Gc122 
                                                 + x*z*cos(2*theta)*kd->Gc122 + 2*x*z*kd->Gc133 - 
                                                 2*x*z*cos(2*theta)*kd->Gc133 + y*z*kd->Gc211 + 
                                                 y*z*cos(2*theta)*kd->Gc211 + y*z*kd->Gc222 + 
                                                 y*z*cos(2*theta)*kd->Gc222 + 2*y*z*kd->Gc233 - 
                                                 2*y*z*cos(2*theta)*kd->Gc233 - pw2(x)*kd->Gc311 - 
                                                 pw2(y)*kd->Gc311 - pw2(x)*cos(2*theta)*kd->Gc311 - 
                                                 pw2(y)*cos(2*theta)*kd->Gc311 + 
                                                 4*pw2(ct)*sin(2*phi)*(x*z*kd->Gc112 + y*z*kd->Gc212 - (pw2(x) 
                                                                                                        + pw2(y))*kd->Gc312) + 
                                                 4*cp*sin(2*theta)*(-(z*(x*kd->Gc113 + y*kd->Gc213)) + 
                                                                    (pw2(x) + pw2(y))*kd->Gc313) - pw2(x)*kd->Gc322 - 
                                                 pw2(y)*kd->Gc322 - pw2(x)*cos(2*theta)*kd->Gc322 - 
                                                 pw2(y)*cos(2*theta)*kd->Gc322 - 
                                                 2*pw2(ct)*cos(2*phi)*(x*z*kd->Gc122 - y*z*kd->Gc211 + 
                                                                       y*z*kd->Gc222 + pw2(x)*kd->Gc311 + pw2(y)*kd->Gc311 - 
                                                                       (pw2(x) + pw2(y))*kd->Gc322) + 
                                                 4*sp*sin(2*theta)*(-(z*(x*kd->Gc123 + y*kd->Gc223)) + 
                                                                    (pw2(x) + pw2(y))*kd->Gc323) - 4*pw2(st)*(pw2(x) + 
                                                                                                              pw2(y))*kd->Gc333))/
          (4.*pow(pw2(r),2.5)*pow((pw2(x) + 
                                   pw2(y))/pw2(r),1.5));
 
        G122[theta_i][phi_i] = (pw2(st)*(z*((pw2(x) + pw2(y))*(pw2(x) + pw2(y) - 
                                                               pw2(z)) - (3*(pw2(x) + pw2(y)) + pw2(z))*((x - y)*(x 
                                                                                                                  + y)*cos(2*phi) + 2*x*y*sin(2*phi))) + 
                                         pw2(r)*(pw2(x) + 
                                                 pw2(y))*(2*pw2(sp)*x*z*kd->Gc111 + x*z*kd->Gc122 + 
                                                          y*z*kd->Gc211 + y*z*kd->Gc222 - pw2(x)*kd->Gc311 - 
                                                          pw2(y)*kd->Gc311 + 
                                                          2*sin(2*phi)*(-(z*(x*kd->Gc112 + y*kd->Gc212)) + 
                                                                        (pw2(x) + pw2(y))*kd->Gc312) - pw2(x)*kd->Gc322 - 
                                                          pw2(y)*kd->Gc322 + 
                                                          cos(2*phi)*(x*z*kd->Gc122 - y*z*kd->Gc211 + y*z*kd->Gc222 + 
                                                                      pw2(x)*kd->Gc311 + pw2(y)*kd->Gc311 - (pw2(x) + 
                                                                                                             pw2(y))*kd->Gc322))))/
          (2.*pow(pw2(r),2.5)*pow((pw2(x) + 
                                   pw2(y))/pw2(r),1.5));
  

        G112[theta_i][phi_i] = (2*(-(sp*x) + cp*y)*(pw2(x) + pw2(y))*(pw2(x) + 
                                                                      pw2(y) - pw2(z)) + 2*(sp*x - cp*y)*(pw2(x) + 
                                                                                                          pw2(y))*(pw2(x) + pw2(y) - pw2(z))*cos(2*theta) + 
                                z*(3*(pw2(x) + pw2(y)) + 
                                   pw2(z))*sin(2*theta)*(2*x*y*cos(2*phi) + (-pw2(x) + 
                                                                             pw2(y))*sin(2*phi)) + 
                                pw2(r)*(pw2(x) + 
                                        pw2(y))*(2*cos(2*phi)*sin(2*theta)*(x*z*kd->Gc112 + y*z*kd->Gc212 
                                                                            - (pw2(x) + pw2(y))*kd->Gc312) + 
                                                 4*sp*pw2(st)*(x*z*kd->Gc113 + y*z*kd->Gc213 - (pw2(x) 
                                                                                                + pw2(y))*kd->Gc313) + 
                                                 sin(2*theta)*sin(2*phi)*(-(x*z*kd->Gc111) + x*z*kd->Gc122 - 
                                                                          y*z*kd->Gc211 + y*z*kd->Gc222 + pw2(x)*kd->Gc311 + 
                                                                          pw2(y)*kd->Gc311 - (pw2(x) + pw2(y))*kd->Gc322) + 
                                                 4*cp*pw2(st)*(-(z*(x*kd->Gc123 + y*kd->Gc223)) + 
                                                               (pw2(x) + 
                                                                pw2(y))*kd->Gc323)))/(4.*pow(pw2(r),2.5)*pow((pw2(x) 
                                                                                                              + pw2(y))/pw2(r),1.5));


        G211[theta_i][phi_i] = (pw2(r)*(4*pw2(ct)*(-2*x*y*cos(2*phi) + (x - y)*(x + 
                                                                                y)*sin(2*phi)) + (pw2(x) + pw2(y))*
                                        (-4*pw2(cp)*pw2(ct)*y*kd->Gc111 - y*kd->Gc122 - 
                                         y*cos(2*theta)*kd->Gc122 - 2*y*kd->Gc133 + 2*y*cos(2*theta)*kd->Gc133 
                                         + x*kd->Gc211 + x*cos(2*theta)*kd->Gc211 + 
                                         4*pw2(ct)*sin(2*phi)*(-(y*kd->Gc112) + x*kd->Gc212) + 
                                         4*cp*sin(2*theta)*(y*kd->Gc113 - x*kd->Gc213) + 
                                         2*pw2(ct)*cos(2*phi)*(y*kd->Gc122 + x*(kd->Gc211 - kd->Gc222)) + 
                                         x*kd->Gc222 + 
                                         x*cos(2*theta)*kd->Gc222 + 4*sp*sin(2*theta)*(y*kd->Gc123 - 
                                                                                       x*kd->Gc223) + 2*x*kd->Gc233 - 
                                         2*x*cos(2*theta)*kd->Gc233)))/(4.*pow(pw2(x) + pw2(y),2));


        G212[theta_i][phi_i] = (pw2(r)*(2*sin(2*theta)*((x - y)*(x + y)*cos(2*phi) + 
                                                        2*x*y*sin(2*phi)) + (pw2(x) + pw2(y))*
                                        (2*cos(2*phi)*sin(2*theta)*(-(y*kd->Gc112) + x*kd->Gc212) + 
                                         4*sp*pw2(st)*(-(y*kd->Gc113) + x*kd->Gc213) + 
                                         sin(2*theta)*sin(2*phi)*(y*kd->Gc111 - y*kd->Gc122 + x*(-kd->Gc211 + 
                                                                                                 kd->Gc222)) + 
                                         4*cp*pw2(st)*(y*kd->Gc123 - 
                                                       x*kd->Gc223))))/(4.*pow(pw2(x) + pw2(y),2));
  

        G222[theta_i][phi_i] = (pw2(r)*pw2(st)*(4*x*y*cos(2*phi) + 2*(-pw2(x) + 
                                                                      pw2(y))*sin(2*phi) - 
                                                (pw2(x) + pw2(y))*(2*pw2(sp)*y*kd->Gc111 - 
                                                                   2*y*sin(2*phi)*kd->Gc112 + y*kd->Gc122 + y*cos(2*phi)*kd->Gc122 - 
                                                                   x*kd->Gc211 + x*cos(2*phi)*kd->Gc211 + 2*x*sin(2*phi)*kd->Gc212 - 
                                                                   2*pw2(cp)*x*kd->Gc222)))/(2.*pow(pw2(x) + 
                                                                                                    pw2(y),2));

        break;
      }
    }  
    if(p != level->end())
      break;

  }

  mpi.Barrier();

  if (mpi.getSize() > 1)
  {
    mpi.AllReduce(&cur_mpi_level, 1, MPI_MAX);
  }
  mpi.Barrier();

  if(mpi.getSize() > 1 && ln != cur_mpi_level)
  {
    cur_mpi_rank = -1;
  }

  mpi.Barrier();

  if (mpi.getSize() > 1)
  {
    mpi.AllReduce(&cur_mpi_rank, 1, MPI_MAX);
  }
  mpi.Barrier();
  
  if(cur_mpi_rank == -1)
    TBOX_ERROR("Cannot find patch cover the cell\n");
  
  mpi.Barrier();
  if (mpi.getSize() > 1 ) {
    mpi.Bcast(&G111[theta_i][phi_i], 1, MPI_DOUBLE, cur_mpi_rank);
    mpi.Bcast(&G112[theta_i][phi_i], 1, MPI_DOUBLE, cur_mpi_rank);
    mpi.Bcast(&G122[theta_i][phi_i], 1, MPI_DOUBLE, cur_mpi_rank);
    mpi.Bcast(&G211[theta_i][phi_i], 1, MPI_DOUBLE, cur_mpi_rank);
    mpi.Bcast(&G212[theta_i][phi_i], 1, MPI_DOUBLE, cur_mpi_rank);
    mpi.Bcast(&G222[theta_i][phi_i], 1, MPI_DOUBLE, cur_mpi_rank);

  }
  mpi.Barrier();

}


void Horizon::set_kd_values(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double theta, double phi, int theta_i, int phi_i, double r, KillingData *kd, BSSN * bssn)
{
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  cur_mpi_rank = -1;
  //if(patch_work_mpi_rank == mpi.getRank())
  int ln_num = hierarchy->getNumberOfLevels(), cur_mpi_level = -1, ln;

  real_t x = r * cos(phi) * sin(theta);
  real_t y = r * sin(phi) * sin(theta);
  real_t z = r * cos(theta);
        
  real_t st = sin(theta);
  real_t ct = cos(theta);
  real_t sp = sin(phi);
  real_t cp = cos(phi);

  for (ln = ln_num - 1; ln >= 0; ln--)
  {
    boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    hier::PatchLevel::iterator p(level->begin());
    for (;p != level->end(); ++p)
    {
      const boost::shared_ptr<hier::Patch>& patch = *p;
      const hier::Box& box = patch->getBox();

      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];
      
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));


      const real_t * dx = &(patch_geometry->getDx())[0];
      int i0 = floor((x + origin[0] - domain_lower[0] ) / dx[0] - 0.5);
      int j0 = floor((y + origin[1] - domain_lower[1] ) / dx[1] - 0.5);
      int k0 = floor((z + origin[2] - domain_lower[2] ) / dx[2] - 0.5);

      
      if( i0 >= lower[0] && i0 <= upper[0]
          && j0 >= lower[1] && j0 <= upper[1]
          && k0 >= lower[2] && k0 <= upper[2])
      {
        cur_mpi_rank = mpi.getRank();

        cur_mpi_level = ln;
        
        initPData(patch);
        initMDA(patch);
        bssn->initPData(patch);
        bssn->initMDA(patch);

        BSSNData bd = {0};


        for(int k = k0; k <= k0 + 1; k++)
        {
          real_t z0 = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0 - origin[2];

          COSMO_APPLY_TO_IJK_PERMS(HORIZON_DEFINE_TEMP_GJ);
          COSMO_APPLY_TO_IJ_PERMS(HORIZON_DEFINE_TEMP_MJ);
          HORIZON_DEFINE_TEMP_RJ;
          for(int j = j0; j <= j0 + 1; j++)
          {
            real_t y0 = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0 - origin[1];

            COSMO_APPLY_TO_IJK_PERMS(HORIZON_DEFINE_TEMP_GI);
            COSMO_APPLY_TO_IJ_PERMS(HORIZON_DEFINE_TEMP_MI);
            HORIZON_DEFINE_TEMP_RI;
            for(int i = i0; i <= i0 + 1; i++)
            {
              real_t x0 = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0 - origin[0];
              bssn->set_bd_values(i, j, k, &bd, dx);                
                          
              COSMO_APPLY_TO_IJK_PERMS(HORIZON_INTERPOLATE_G_1);
              COSMO_APPLY_TO_IJ_PERMS(HORIZON_INTERPOLATE_M_1);
              HORIZON_INTERPOLATE_R_1;
            }
      
            COSMO_APPLY_TO_IJK_PERMS(HORIZON_INTERPOLATE_G_2);
            COSMO_APPLY_TO_IJ_PERMS(HORIZON_INTERPOLATE_M_2);
            HORIZON_INTERPOLATE_R_2;
          }
    
          COSMO_APPLY_TO_IJK_PERMS(HORIZON_INTERPOLATE_G_3);
          COSMO_APPLY_TO_IJ_PERMS(HORIZON_INTERPOLATE_M_3);
          HORIZON_INTERPOLATE_R_3;
        }

  
        kd->q11 = pw2(r)*(ct*(pw2(cp)*ct*kd->m11
                              + ct*sin(2*phi)*kd->m12
                              + ct*pw2(sp)*kd->m22
                              - 2*st*(cp*kd->m13 + sp*kd->m23)) + pw2(st)*kd->m33);

  
        kd->q12 = pw2(r)*st*(ct*cos(2*phi)*kd->m12
                             + sp*st*kd->m13
                             + cp*ct*sp*(-kd->m11 + kd->m22) - cp*st*kd->m23);

        kd->q22 = pw2(r)*pw2(st)*(pw2(sp)*kd->m11 + cp*(-2*sp*kd->m12 + cp*kd->m22));

        kd->Gs111 = ((pw2(x) + pw2(y))*z*(5*(pw2(x) + pw2(y)) - 
                                          pw2(z))*cos(2*theta) + 
                     (pw2(x) + pw2(y))*(-(z*(3*(pw2(x) + pw2(y)) + 
                                             pw2(z))) + 4*(cp*x + sp*y)*(pw2(x) + pw2(y) - 
                                                                         pw2(z))*sin(2*theta)) + 
                     2*pw2(ct)*z*(3*(pw2(x) + pw2(y)) + pw2(z))*((x - 
                                                                  y)*(x + y)*cos(2*phi) + 2*x*y*sin(2*phi)) + 
                     pw2(r)*(pw2(x) + 
                             pw2(y))*(4*pw2(cp)*pw2(ct)*x*z*kd->Gc111 + x*z*kd->Gc122 
                                      + x*z*cos(2*theta)*kd->Gc122 + 2*x*z*kd->Gc133 - 
                                      2*x*z*cos(2*theta)*kd->Gc133 + y*z*kd->Gc211 + 
                                      y*z*cos(2*theta)*kd->Gc211 + y*z*kd->Gc222 + 
                                      y*z*cos(2*theta)*kd->Gc222 + 2*y*z*kd->Gc233 - 
                                      2*y*z*cos(2*theta)*kd->Gc233 - pw2(x)*kd->Gc311 - 
                                      pw2(y)*kd->Gc311 - pw2(x)*cos(2*theta)*kd->Gc311 - 
                                      pw2(y)*cos(2*theta)*kd->Gc311 + 
                                      4*pw2(ct)*sin(2*phi)*(x*z*kd->Gc112 + y*z*kd->Gc212 - (pw2(x) 
                                                                                             + pw2(y))*kd->Gc312) + 
                                      4*cp*sin(2*theta)*(-(z*(x*kd->Gc113 + y*kd->Gc213)) + 
                                                         (pw2(x) + pw2(y))*kd->Gc313) - pw2(x)*kd->Gc322 - 
                                      pw2(y)*kd->Gc322 - pw2(x)*cos(2*theta)*kd->Gc322 - 
                                      pw2(y)*cos(2*theta)*kd->Gc322 - 
                                      2*pw2(ct)*cos(2*phi)*(x*z*kd->Gc122 - y*z*kd->Gc211 + 
                                                            y*z*kd->Gc222 + pw2(x)*kd->Gc311 + pw2(y)*kd->Gc311 - 
                                                            (pw2(x) + pw2(y))*kd->Gc322) + 
                                      4*sp*sin(2*theta)*(-(z*(x*kd->Gc123 + y*kd->Gc223)) + 
                                                         (pw2(x) + pw2(y))*kd->Gc323) - 4*pw2(st)*(pw2(x) + 
                                                                                                   pw2(y))*kd->Gc333))/
          (4.*pow(pw2(r),2.5)*pow((pw2(x) + 
                                   pw2(y))/pw2(r),1.5));
 
        kd->Gs122 = (pw2(st)*(z*((pw2(x) + pw2(y))*(pw2(x) + pw2(y) - 
                                                    pw2(z)) - (3*(pw2(x) + pw2(y)) + pw2(z))*((x - y)*(x 
                                                                                                       + y)*cos(2*phi) + 2*x*y*sin(2*phi))) + 
                              pw2(r)*(pw2(x) + 
                                      pw2(y))*(2*pw2(sp)*x*z*kd->Gc111 + x*z*kd->Gc122 + 
                                               y*z*kd->Gc211 + y*z*kd->Gc222 - pw2(x)*kd->Gc311 - 
                                               pw2(y)*kd->Gc311 + 
                                               2*sin(2*phi)*(-(z*(x*kd->Gc112 + y*kd->Gc212)) + 
                                                             (pw2(x) + pw2(y))*kd->Gc312) - pw2(x)*kd->Gc322 - 
                                               pw2(y)*kd->Gc322 + 
                                               cos(2*phi)*(x*z*kd->Gc122 - y*z*kd->Gc211 + y*z*kd->Gc222 + 
                                                           pw2(x)*kd->Gc311 + pw2(y)*kd->Gc311 - (pw2(x) + 
                                                                                                  pw2(y))*kd->Gc322))))/
          (2.*pow(pw2(r),2.5)*pow((pw2(x) + 
                                   pw2(y))/pw2(r),1.5));
  

        kd->Gs112 = (2*(-(sp*x) + cp*y)*(pw2(x) + pw2(y))*(pw2(x) + 
                                                           pw2(y) - pw2(z)) + 2*(sp*x - cp*y)*(pw2(x) + 
                                                                                               pw2(y))*(pw2(x) + pw2(y) - pw2(z))*cos(2*theta) + 
                     z*(3*(pw2(x) + pw2(y)) + 
                        pw2(z))*sin(2*theta)*(2*x*y*cos(2*phi) + (-pw2(x) + 
                                                                  pw2(y))*sin(2*phi)) + 
                     pw2(r)*(pw2(x) + 
                             pw2(y))*(2*cos(2*phi)*sin(2*theta)*(x*z*kd->Gc112 + y*z*kd->Gc212 
                                                                 - (pw2(x) + pw2(y))*kd->Gc312) + 
                                      4*sp*pw2(st)*(x*z*kd->Gc113 + y*z*kd->Gc213 - (pw2(x) 
                                                                                     + pw2(y))*kd->Gc313) + 
                                      sin(2*theta)*sin(2*phi)*(-(x*z*kd->Gc111) + x*z*kd->Gc122 - 
                                                               y*z*kd->Gc211 + y*z*kd->Gc222 + pw2(x)*kd->Gc311 + 
                                                               pw2(y)*kd->Gc311 - (pw2(x) + pw2(y))*kd->Gc322) + 
                                      4*cp*pw2(st)*(-(z*(x*kd->Gc123 + y*kd->Gc223)) + 
                                                    (pw2(x) + 
                                                     pw2(y))*kd->Gc323)))/(4.*pow(pw2(r),2.5)*pow((pw2(x) 
                                                                                                   + pw2(y))/pw2(r),1.5));


        kd->Gs211 = (pw2(r)*(4*pw2(ct)*(-2*x*y*cos(2*phi) + (x - y)*(x + 
                                                                     y)*sin(2*phi)) + (pw2(x) + pw2(y))*
                             (-4*pw2(cp)*pw2(ct)*y*kd->Gc111 - y*kd->Gc122 - 
                              y*cos(2*theta)*kd->Gc122 - 2*y*kd->Gc133 + 2*y*cos(2*theta)*kd->Gc133 
                              + x*kd->Gc211 + x*cos(2*theta)*kd->Gc211 + 
                              4*pw2(ct)*sin(2*phi)*(-(y*kd->Gc112) + x*kd->Gc212) + 
                              4*cp*sin(2*theta)*(y*kd->Gc113 - x*kd->Gc213) + 
                              2*pw2(ct)*cos(2*phi)*(y*kd->Gc122 + x*(kd->Gc211 - kd->Gc222)) + 
                              x*kd->Gc222 + 
                              x*cos(2*theta)*kd->Gc222 + 4*sp*sin(2*theta)*(y*kd->Gc123 - 
                                                                            x*kd->Gc223) + 2*x*kd->Gc233 - 
                              2*x*cos(2*theta)*kd->Gc233)))/(4.*pow(pw2(x) + pw2(y),2));


        kd->Gs212 = (pw2(r)*(2*sin(2*theta)*((x - y)*(x + y)*cos(2*phi) + 
                                             2*x*y*sin(2*phi)) + (pw2(x) + pw2(y))*
                             (2*cos(2*phi)*sin(2*theta)*(-(y*kd->Gc112) + x*kd->Gc212) + 
                              4*sp*pw2(st)*(-(y*kd->Gc113) + x*kd->Gc213) + 
                              sin(2*theta)*sin(2*phi)*(y*kd->Gc111 - y*kd->Gc122 + x*(-kd->Gc211 + 
                                                                                      kd->Gc222)) + 
                              4*cp*pw2(st)*(y*kd->Gc123 - 
                                            x*kd->Gc223))))/(4.*pow(pw2(x) + pw2(y),2));
  

        kd->Gs222 = (pw2(r)*pw2(st)*(4*x*y*cos(2*phi) + 2*(-pw2(x) + 
                                                           pw2(y))*sin(2*phi) - 
                                     (pw2(x) + pw2(y))*(2*pw2(sp)*y*kd->Gc111 - 
                                                        2*y*sin(2*phi)*kd->Gc112 + y*kd->Gc122 + y*cos(2*phi)*kd->Gc122 - 
                                                        x*kd->Gc211 + x*cos(2*phi)*kd->Gc211 + 2*x*sin(2*phi)*kd->Gc212 - 
                                                        2*pw2(cp)*x*kd->Gc222)))/(2.*pow(pw2(x) + 
                                                                                         pw2(y),2));
  
        double det = kd->q11 * kd->q22 - kd->q12 * kd->q12;


        kd->qi11 = kd->q22 / det;
        kd->qi12 = -kd->q12 / det;
        kd->qi22 =  kd->q11 / det;


        kd->R11 = HORIZON_CALCULATE_D1G(1,1,1)
          + HORIZON_CALCULATE_D2G(2,1,1)
          - HORIZON_CALCULATE_D1G(1,1,1)
          - HORIZON_CALCULATE_D1G(2,2,1)
          + (kd->Gs111 * kd->Gs111 + kd->Gs221 * kd->Gs111
             + kd->Gs112 * kd->Gs211 + kd->Gs222 * kd->Gs211)
          - (kd->Gs111 * kd->Gs111 + kd->Gs112 * kd->Gs211
             + kd->Gs211 * kd->Gs121 + kd->Gs212 * kd->Gs221);
  
        kd->R12 = HORIZON_CALCULATE_D1G(1,1,2)
          + HORIZON_CALCULATE_D2G(2,1,2)
          - HORIZON_CALCULATE_D2G(1,1,1)
          - HORIZON_CALCULATE_D2G(2,2,1)
          + (kd->Gs111 * kd->Gs112 + kd->Gs221 * kd->Gs112
             + kd->Gs112 * kd->Gs212 + kd->Gs222 * kd->Gs212)
          - (kd->Gs121 * kd->Gs121 + kd->Gs122 * kd->Gs211
             + kd->Gs221 * kd->Gs121 + kd->Gs222 * kd->Gs221);

        kd->R22 = HORIZON_CALCULATE_D1G(1,2,2)
          + HORIZON_CALCULATE_D2G(2,2,2)
          - HORIZON_CALCULATE_D2G(1,1,2)
          - HORIZON_CALCULATE_D2G(2,2,2)
          + (kd->Gs111 * kd->Gs122 + kd->Gs221 * kd->Gs122
             + kd->Gs112 * kd->Gs222 + kd->Gs222 * kd->Gs222)
          - (kd->Gs121 * kd->Gs122 + kd->Gs122 * kd->Gs212
             + kd->Gs221 * kd->Gs122 + kd->Gs222 * kd->Gs222);

        kd->R = kd->qi11 * kd->R11 + kd->qi12 * kd->R12 * 2.0
          + kd->qi22 * kd->R22;

        break;
      }
    }
    if(p != level->end())
      break;

  }

  mpi.Barrier();

  if (mpi.getSize() > 1)
  {
    mpi.AllReduce(&cur_mpi_level, 1, MPI_MAX);
  }
  mpi.Barrier();

  if(mpi.getSize() > 1 && ln != cur_mpi_level)
  {
    cur_mpi_rank = -1;
  }

  
  mpi.Barrier();

  if (mpi.getSize() > 1)
  {
    mpi.AllReduce(&cur_mpi_rank, 1, MPI_MAX);
  }
  mpi.Barrier();
  
  if(cur_mpi_rank == -1)
    TBOX_ERROR("Cannot find proper patch in set bd values\n");

}

real_t Horizon::ev_k_theta_dphi(KillingData *kd, int theta_i, int phi_i)
{
  return kd->Gs112 * k_theta[theta_i][phi_i] + kd->Gs212 * k_phi[theta_i][phi_i]
    - k_L[theta_i][phi_i] * sqrt(kd->q11 * kd->q22 - kd->q12 * kd->q12);
}

real_t Horizon::ev_k_phi_dphi(KillingData *kd, int theta_i, int phi_i)
{
  return kd->Gs122 * k_theta[theta_i][phi_i] + kd->Gs222 * k_phi[theta_i][phi_i];
}

real_t Horizon::ev_k_L_dphi(KillingData *kd, int theta_i, int phi_i)
{
  return 0.5 * kd->R * sqrt(kd->q11 * kd->q22 - kd->q12 * kd->q12)
    * (kd->qi11 * k_theta[theta_i][phi_i] + kd->qi12 * k_phi[theta_i][phi_i]);
}

real_t Horizon::ev_k_theta_dtheta(KillingData *kd, int theta_i, int phi_i)
{
  return kd->Gs111 * k_theta[theta_i][phi_i] + kd->Gs211 * k_phi[theta_i][phi_i];
}

real_t Horizon::ev_k_phi_dtheta(KillingData *kd, int theta_i, int phi_i)
{
  return kd->Gs112 * k_theta[theta_i][phi_i] + kd->Gs212 * k_phi[theta_i][phi_i]
    + k_L[theta_i][phi_i] * sqrt(kd->q11 * kd->q22 - kd->q12 * kd->q12);
}

real_t Horizon::ev_k_L_dtheta(KillingData *kd, int theta_i, int phi_i)
{
  return -0.5 * kd->R * sqrt(kd->q11 * kd->q22 - kd->q12 * kd->q12)
        * (kd->qi12 * k_theta[theta_i][phi_i] + kd->qi22 * k_phi[theta_i][phi_i]);
    //* (kd->qi11 * k_theta[theta_i][phi_i] + kd->qi12 * k_phi[theta_i][phi_i]);
}

real_t Horizon::getRadius(int theta_i, int phi_i)
{
  double res = 0;
  if(theta_i + 1 >= 2 * n_theta || theta_i < 0)
    TBOX_ERROR("Theta_i is out of the range!\b");
  
  for(int i = theta_i; i <= theta_i + 1 ; i++)
  {
    for(int j = phi_i; j <= phi_i + 1; j ++)
    {
      res += ah_radius[i][(j+2*n_phi)%(2*n_phi)];
    }
  }
  return res / 4.0;
}
  
void Horizon::transportKillingTheta(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t phi_i, double k_theta_0, double k_phi_0, double k_L_0, BSSN * bssn)
{
  // transporting killing vector (or test vector from 0 to 2 \pi)
  double phi = 2.0 * PI * ((double) phi_i +0.5) / (double)n_phi;
  double dtheta = PI / (double) n_theta;

  // always starts from n_theta / 2
  k_theta[n_theta/2][phi_i] = k_theta_0;
  k_phi[n_theta/2][phi_i] = k_phi_0;
  k_L[n_theta/2][phi_i] = k_L_0;

  KillingData kd = {0};

  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());

  int i_bak = patch_work_i;
  int j_bak = patch_work_j;
  int k_bak = patch_work_k;
  int mpi_bak = patch_work_mpi_rank;
  int local_id_bak = local_id;
  int level_bak = patch_work_level;
  
  // transporting theta_i to theta_i + 1
  for(int theta_i = n_theta/2; theta_i < n_theta - 1; theta_i++)
  {
    real_t k_phi_0 = k_phi[theta_i][phi_i];
    real_t k_theta_0 = k_theta[theta_i][phi_i];
    real_t k_L_0 = k_L[theta_i][phi_i];
    kd = {0};
    
    /********Doing K1 ******************************************/
    double theta =  PI * ((double) theta_i +0.5) / (double)n_theta;
    mpi.Barrier();
    //    real_t r = ah_radius[theta_i*2][phi_i*2];
    real_t r = getRadius(theta_i*2, phi_i*2);
    mpi.Barrier();
    
    set_kd_values(hierarchy, theta, phi, theta_i*2, phi_i*2, r, &kd, bssn);
    
    real_t k1_theta = ev_k_theta_dtheta(&kd, theta_i, phi_i);
    real_t k1_phi = ev_k_phi_dtheta(&kd, theta_i, phi_i);
    real_t k1_L = ev_k_L_dtheta(&kd, theta_i, phi_i);

    mpi.Barrier();
    if (mpi.getSize() > 1 ) {
      mpi.Bcast(&k1_theta, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k1_phi, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k1_L, 1, MPI_DOUBLE, cur_mpi_rank);
    }
    mpi.Barrier();
    kd = {0};
    /********Doing K2 ******************************************/
    theta += dtheta / 2.0;

    k_theta[theta_i][phi_i] = k_theta_0 + dtheta * k1_theta / 2.0;
    k_phi[theta_i][phi_i] = k_phi_0 + dtheta * k1_phi / 2.0;    
    k_L[theta_i][phi_i] = k_L_0 + dtheta * k1_L / 2.0;


    mpi.Barrier();
    //    r =  findRadius(hierarchy, theta, phi);
    //    r = ah_radius[theta_i * 2 + 1][phi_i*2];
    r = getRadius(theta_i * 2 +1, phi_i *2);
    
    mpi.Barrier();
    
    set_kd_values(hierarchy, theta, phi, (theta_i*2 +1)%(2*n_theta), phi_i*2, r, &kd, bssn);

    real_t k2_theta = ev_k_theta_dtheta(&kd, theta_i, phi_i);
    real_t k2_phi = ev_k_phi_dtheta(&kd, theta_i, phi_i);
    real_t k2_L = ev_k_L_dtheta(&kd, theta_i, phi_i);

    mpi.Barrier();
    if (mpi.getSize() > 1 ) {
      mpi.Bcast(&k2_theta, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k2_phi, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k2_L, 1, MPI_DOUBLE, cur_mpi_rank);
    }
    mpi.Barrier();
    /********Doing K3 ******************************************/
    k_theta[theta_i][phi_i] = k_theta_0 + dtheta * k2_theta / 2.0;
    k_phi[theta_i][phi_i] = k_phi_0 + dtheta * k2_phi / 2.0;    
    k_L[theta_i][phi_i] = k_L_0 + dtheta * k2_L / 2.0;


    
    real_t k3_theta = ev_k_theta_dtheta(&kd, theta_i, phi_i);
    real_t k3_phi = ev_k_phi_dtheta(&kd, theta_i, phi_i);
    real_t k3_L = ev_k_L_dtheta(&kd, theta_i, phi_i);

    mpi.Barrier();
    if (mpi.getSize() > 1 ) {
      mpi.Bcast(&k3_theta, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k3_phi, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k3_L, 1, MPI_DOUBLE, cur_mpi_rank);
    }
    mpi.Barrier();
    kd = {0};
    /********Doing K4 ******************************************/

    theta += dtheta / 2.0;

    k_theta[theta_i][phi_i] = k_theta_0 + dtheta * k3_theta;
    k_phi[theta_i][phi_i] = k_phi_0 + dtheta * k3_phi;    
    k_L[theta_i][phi_i] = k_L_0 + dtheta * k3_L;
    mpi.Barrier();
    //    r =  findRadius(hierarchy, theta, phi);
    //    r = ah_radius[theta_i * 2 + 2][phi_i*2];
    r = getRadius(theta_i * 2 +2, phi_i *2);
    mpi.Barrier();
    set_kd_values(hierarchy, theta, phi, (theta_i*2 +2)%(2*n_theta), phi_i*2, r, &kd, bssn);


    
    real_t k4_theta = ev_k_theta_dtheta(&kd, theta_i, phi_i);
    real_t k4_phi = ev_k_phi_dtheta(&kd, theta_i, phi_i);
    real_t k4_L = ev_k_L_dtheta(&kd, theta_i, phi_i);

    mpi.Barrier();
    if (mpi.getSize() > 1 ) {
      mpi.Bcast(&k4_theta, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k4_phi, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k4_L, 1, MPI_DOUBLE, cur_mpi_rank);
    }
    mpi.Barrier();
    
    k_theta[(theta_i+1)%n_theta][phi_i] = k_theta_0 + dtheta /6.0 *
      ( k1_theta + 2.0 * k2_theta + 2.0 * k3_theta + k4_theta);
    k_phi[(theta_i+1)%n_theta][phi_i] = k_phi_0 + dtheta /6.0 *
      ( k1_phi + 2.0 * k2_phi + 2.0 * k3_phi + k4_phi);
    k_L[(theta_i+1)%n_theta][phi_i] = k_L_0 + dtheta /6.0 *
      ( k1_L + 2.0 * k2_L + 2.0 * k3_L + k4_L);

    // restore the initial result
    k_phi[theta_i][phi_i] = k_phi_0;
    k_theta[theta_i][phi_i] = k_theta_0;
    k_L[theta_i][phi_i] = k_L_0;

  }

  // restore the original choice for i, j, k
  // and evolve killing field downward along theta
  patch_work_i = i_bak;
  patch_work_j = j_bak;
  patch_work_k = k_bak;
  patch_work_mpi_rank = mpi_bak;
  local_id = local_id_bak;
  patch_work_level = level_bak;

  
  
  dtheta = - dtheta;
  // transporting theta_i to theta_i - 1
  for(int theta_i = n_theta/2; theta_i > 0; theta_i--)
  {
    real_t k_phi_0 = k_phi[theta_i][phi_i];
    real_t k_theta_0 = k_theta[theta_i][phi_i];
    real_t k_L_0 = k_L[theta_i][phi_i];
    kd = {0};
    /********Doing K1 ******************************************/
    double theta =  PI * ((double) theta_i +0.5) / (double)n_theta;
    mpi.Barrier();
    //    real_t r = ah_radius[theta_i*2][phi_i*2];
    real_t r = getRadius(theta_i*2, phi_i*2);

    mpi.Barrier();
    set_kd_values(hierarchy, theta, phi, theta_i*2, phi_i*2, r, &kd, bssn);
    
    real_t k1_theta = ev_k_theta_dtheta(&kd, theta_i, phi_i);
    real_t k1_phi = ev_k_phi_dtheta(&kd, theta_i, phi_i);
    real_t k1_L = ev_k_L_dtheta(&kd, theta_i, phi_i);

    mpi.Barrier();
    if (mpi.getSize() > 1 ) {
      mpi.Bcast(&k1_theta, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k1_phi, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k1_L, 1, MPI_DOUBLE, cur_mpi_rank);
    }
    mpi.Barrier();

    kd = {0};
    /********Doing K2 ******************************************/
    theta += dtheta / 2.0;

    k_theta[theta_i][phi_i] = k_theta_0 + dtheta * k1_theta / 2.0;
    k_phi[theta_i][phi_i] = k_phi_0 + dtheta * k1_phi / 2.0;    
    k_L[theta_i][phi_i] = k_L_0 + dtheta * k1_L / 2.0;


    mpi.Barrier();
    //    r =  findRadius(hierarchy, theta, phi);
    //    r = ah_radius[theta_i*2-1][phi_i*2];
    r = getRadius(theta_i*2-1, phi_i*2);
    mpi.Barrier();
    set_kd_values(hierarchy, theta, phi, (theta_i*2 +1)%(2*n_theta), phi_i*2, r, &kd, bssn);

    real_t k2_theta = ev_k_theta_dtheta(&kd, theta_i, phi_i);
    real_t k2_phi = ev_k_phi_dtheta(&kd, theta_i, phi_i);
    real_t k2_L = ev_k_L_dtheta(&kd, theta_i, phi_i);

    mpi.Barrier();
    if (mpi.getSize() > 1 ) {
      mpi.Bcast(&k2_theta, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k2_phi, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k2_L, 1, MPI_DOUBLE, cur_mpi_rank);
    }
    mpi.Barrier();
        
    /********Doing K3 ******************************************/
    k_theta[theta_i][phi_i] = k_theta_0 + dtheta * k2_theta / 2.0;
    k_phi[theta_i][phi_i] = k_phi_0 + dtheta * k2_phi / 2.0;    
    k_L[theta_i][phi_i] = k_L_0 + dtheta * k2_L / 2.0;


    
    real_t k3_theta = ev_k_theta_dtheta(&kd, theta_i, phi_i);
    real_t k3_phi = ev_k_phi_dtheta(&kd, theta_i, phi_i);
    real_t k3_L = ev_k_L_dtheta(&kd, theta_i, phi_i);

    mpi.Barrier();
    if (mpi.getSize() > 1 ) {
      mpi.Bcast(&k3_theta, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k3_phi, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k3_L, 1, MPI_DOUBLE, cur_mpi_rank);
    }
    mpi.Barrier();

    kd = {0};
    /********Doing K4 ******************************************/

    theta += dtheta / 2.0;

    k_theta[theta_i][phi_i] = k_theta_0 + dtheta * k3_theta;
    k_phi[theta_i][phi_i] = k_phi_0 + dtheta * k3_phi;    
    k_L[theta_i][phi_i] = k_L_0 + dtheta * k3_L;
    mpi.Barrier();
    //    r =  findRadius(hierarchy, theta, phi);
    //    r = ah_radius[theta_i*2-2][phi_i*2];
    r= getRadius(theta_i*2-2, phi_i*2);
    mpi.Barrier();
    set_kd_values(hierarchy, theta, phi, (theta_i*2 +2)%(2*n_theta), phi_i*2, r, &kd, bssn);


    
    real_t k4_theta = ev_k_theta_dtheta(&kd, theta_i, phi_i);
    real_t k4_phi = ev_k_phi_dtheta(&kd, theta_i, phi_i);
    real_t k4_L = ev_k_L_dtheta(&kd, theta_i, phi_i);

    mpi.Barrier();
    if (mpi.getSize() > 1 ) {
      mpi.Bcast(&k4_theta, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k4_phi, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k4_L, 1, MPI_DOUBLE, cur_mpi_rank);
    }
    mpi.Barrier();
    
    k_theta[(theta_i-1+n_theta)%n_theta][phi_i] = k_theta_0 + dtheta /6.0 *
      ( k1_theta + 2.0 * k2_theta + 2.0 * k3_theta + k4_theta);
    k_phi[(theta_i-1+n_theta)%n_theta][phi_i] = k_phi_0 + dtheta /6.0 *
      ( k1_phi + 2.0 * k2_phi + 2.0 * k3_phi + k4_phi);
    k_L[(theta_i-1+n_theta)%n_theta][phi_i] = k_L_0 + dtheta /6.0 *
      ( k1_L + 2.0 * k2_L + 2.0 * k3_L + k4_L);

    // restore the initial result
    k_phi[theta_i][phi_i] = k_phi_0;
    k_theta[theta_i][phi_i] = k_theta_0;
    k_L[theta_i][phi_i] = k_L_0;

  }

  // restore the original choice for i, j, k
  // and evolve killing field downward along theta
  patch_work_i = i_bak;
  patch_work_j = j_bak;
  patch_work_k = k_bak;
  patch_work_mpi_rank = mpi_bak;
  local_id = local_id_bak;
  patch_work_level = level_bak;
  
}

void Horizon::transportKillingPhi(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t theta_i, idx_t phi_f, double k_theta_0, double k_phi_0, double k_L_0, BSSN * bssn)
{
  // transporting killing vector (or test vector from 0 to 2 \pi)
  double theta = PI * ((double) theta_i +0.5) / (double)n_theta;
  double dphi = 2.0 * PI / (double) n_phi;

  k_theta[theta_i][0] = k_theta_0;
  k_phi[theta_i][0] = k_phi_0;
  k_L[theta_i][0] = k_L_0;

  KillingData kd = {0};

  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());

  // transporting phi_i to phi_f
  for(int phi_i = 0; phi_i < phi_f; phi_i++)
  {
    real_t k_phi_0 = k_phi[theta_i][phi_i];
    real_t k_theta_0 = k_theta[theta_i][phi_i];
    real_t k_L_0 = k_L[theta_i][phi_i];
    kd = {0};
    
    /********Doing K1 ******************************************/
    double phi = 2.0 * PI * ((double) phi_i +0.5) / (double)n_phi;

    mpi.Barrier();
    //    real_t r = ah_radius[theta_i*2][phi_i*2];
    real_t r = getRadius(theta_i*2, phi_i*2);
    mpi.Barrier();
    
    set_kd_values(hierarchy, theta, phi, theta_i*2 , phi_i*2, r, &kd, bssn);


    real_t k1_theta = ev_k_theta_dphi(&kd, theta_i, phi_i);
    real_t k1_phi = ev_k_phi_dphi(&kd, theta_i, phi_i);
    real_t k1_L = ev_k_L_dphi(&kd, theta_i, phi_i);

    mpi.Barrier();
        
    if (mpi.getSize() > 1) {
      mpi.Bcast(&k1_theta, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k1_phi, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k1_L, 1, MPI_DOUBLE, cur_mpi_rank);
    }
    mpi.Barrier();

    kd = {0};

    /********Doing K2 ******************************************/
    phi += dphi / 2.0;

    k_theta[theta_i][phi_i] = k_theta_0 + dphi * k1_theta / 2.0;
    k_phi[theta_i][phi_i] = k_phi_0 + dphi * k1_phi / 2.0;    
    k_L[theta_i][phi_i] = k_L_0 + dphi * k1_L / 2.0;

    mpi.Barrier();


    //    r =  findRadius(hierarchy, theta, phi);
    //    r = ah_radius[theta_i*2][(phi_i*2 + 1)%n_phi];
    r = getRadius(theta_i*2, phi_i*2+1);

    mpi.Barrier();
    set_kd_values(hierarchy, theta, phi, theta_i*2, (phi_i*2+1)%(2*n_phi), r, &kd, bssn);

    real_t k2_theta = ev_k_theta_dphi(&kd, theta_i, phi_i);
    real_t k2_phi = ev_k_phi_dphi(&kd, theta_i, phi_i);
    real_t k2_L = ev_k_L_dphi(&kd, theta_i, phi_i);

    mpi.Barrier();
    if (mpi.getSize() > 1) {
      mpi.Bcast(&k2_theta, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k2_phi, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k2_L, 1, MPI_DOUBLE, cur_mpi_rank);
    }
    mpi.Barrier();


    //kd = {0};    
    /********Doing K3 ******************************************/
    k_theta[theta_i][phi_i] = k_theta_0 + dphi * k2_theta / 2.0;
    k_phi[theta_i][phi_i] = k_phi_0 + dphi * k2_phi / 2.0;    
    k_L[theta_i][phi_i] = k_L_0 + dphi * k2_L / 2.0;


    real_t k3_theta = ev_k_theta_dphi(&kd, theta_i, phi_i);
    real_t k3_phi = ev_k_phi_dphi(&kd, theta_i, phi_i);
    real_t k3_L = ev_k_L_dphi(&kd, theta_i, phi_i);


    mpi.Barrier();
    if (mpi.getSize() > 1 ) {
      mpi.Bcast(&k3_theta, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k3_phi, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k3_L, 1, MPI_DOUBLE, cur_mpi_rank);
    }
    mpi.Barrier();
    kd = {0};
    /********Doing K4 ******************************************/

    phi += dphi / 2.0;

    k_theta[theta_i][phi_i] = k_theta_0 + dphi * k3_theta;
    k_phi[theta_i][phi_i] = k_phi_0 + dphi * k3_phi;    
    k_L[theta_i][phi_i] = k_L_0 + dphi * k3_L;

    mpi.Barrier();
    //    r = ah_radius[theta_i*2][(phi_i*2 + 2)%n_phi];
    r =  getRadius(theta_i*2, phi_i*2+2);


    mpi.Barrier();

    set_kd_values(hierarchy, theta, phi, theta_i*2, (phi_i*2+2)%(2*n_phi), r, &kd, bssn);
    
    real_t k4_theta = ev_k_theta_dphi(&kd, theta_i, phi_i);
    real_t k4_phi = ev_k_phi_dphi(&kd, theta_i, phi_i);
    real_t k4_L = ev_k_L_dphi(&kd, theta_i, phi_i);

    mpi.Barrier();
    if (mpi.getSize() > 1 ) {
      mpi.Bcast(&k4_theta, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k4_phi, 1, MPI_DOUBLE, cur_mpi_rank);
      mpi.Bcast(&k4_L, 1, MPI_DOUBLE, cur_mpi_rank);
    }
    mpi.Barrier();
    
    k_theta[theta_i][(phi_i+1)%n_phi] = k_theta_0 + dphi /6.0 *
      ( k1_theta + 2.0 * k2_theta + 2.0 * k3_theta + k4_theta);
    k_phi[theta_i][(phi_i+1)%n_phi] = k_phi_0 + dphi /6.0 *
      ( k1_phi + 2.0 * k2_phi + 2.0 * k3_phi + k4_phi);
    k_L[theta_i][(phi_i+1)%n_phi] = k_L_0 + dphi /6.0 *
      ( k1_L + 2.0 * k2_L + 2.0 * k3_L + k4_L);


    // restore the initial result
    k_phi[theta_i][phi_i] = k_phi_0;
    k_theta[theta_i][phi_i] = k_theta_0;
    k_L[theta_i][phi_i] = k_L_0;
  }
}

void Horizon::initG(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn)
{

  KillingData kd = {0};

  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());

  double dphi = 2.0 * PI / (double) n_phi / 2;
  double dtheta = PI / (double) n_theta / 2;

  double min_d0 = min_d;
  
  int theta_i = n_theta, phi_i = 0;

  double phi, theta;
  for(phi_i = 0; phi_i < n_phi*2; phi_i++)
  {
    int patch_work_i_bak = patch_work_i,
      patch_work_j_bak = patch_work_j,
      patch_work_k_bak = patch_work_k,
      patch_work_level_bak = patch_work_level,
      local_id_bak = local_id,
      patch_work_mpi_rank_bak = patch_work_mpi_rank;

    phi = 2.0 * PI * ((double) phi_i +0.5) / (double)n_phi / 2.0;
    theta =  PI * ((double) n_theta + 0.5) / (double)n_theta / 2.0;

    theta_i = n_theta;


    mpi.Barrier();
    real_t r =  findRadius(hierarchy, theta, phi);
    mpi.Barrier();

    kd = {0};
    set_G_values(hierarchy, theta, phi, theta_i, phi_i, r, &kd, bssn);
    ah_radius[theta_i][phi_i] = r;
    for(theta_i = n_theta+1; theta_i < 2 * n_theta; theta_i++)
    {
      theta =  PI * ((double) theta_i + 0.5) / (double)n_theta / 2.0;
      mpi.Barrier();

      real_t r =  findRadius(hierarchy, theta, phi);
      mpi.Barrier();

      kd = {0};

      
      set_G_values(hierarchy, theta, phi, theta_i, phi_i, r, &kd, bssn);
      mpi.Barrier();
      ah_radius[theta_i][phi_i] = r;
    }

    patch_work_i = patch_work_i_bak, patch_work_j = patch_work_j_bak, patch_work_k = patch_work_k_bak;
    patch_work_level = patch_work_level_bak, local_id = local_id_bak, patch_work_mpi_rank = patch_work_mpi_rank_bak;
    
    for(theta_i = n_theta - 1; theta_i >= 0; theta_i--)
    {
      theta =  PI * ((double) theta_i + 0.5) / (double)n_theta / 2.0;
      mpi.Barrier();

      real_t r =  findRadius(hierarchy, theta, phi);
      mpi.Barrier();
      kd = {0};
      set_G_values(hierarchy, theta, phi, theta_i, phi_i, r, &kd, bssn);
      mpi.Barrier();
      ah_radius[theta_i][phi_i] = r;
    }
    patch_work_i = patch_work_i_bak, patch_work_j = patch_work_j_bak, patch_work_k = patch_work_k_bak;
    patch_work_level = patch_work_level_bak, local_id = local_id_bak;
    patch_work_mpi_rank = patch_work_mpi_rank_bak;

  }
  phi = 2.0 * PI * ((double) 0.5) / (double)n_phi / 2.0;
  theta =  PI * ((double) n_theta + 0.5) / (double)n_theta / 2.0;
  mpi.Barrier();
  real_t r =  findRadius(hierarchy, theta, phi);
  mpi.Barrier();

  phi = 2.0 * PI * ((double) 0.5) / (double)n_phi;
  theta =  PI * ((double) n_theta/2 + 0.5) / (double)n_theta;
  mpi.Barrier();
  r =  findRadius(hierarchy, PI/2, 0);
  mpi.Barrier();

  if(min_d > min_d0 + EPS)
  {

    TBOX_ERROR("Pointer does not return to the starting patch!\n");

  }


}


void Horizon::findPatch(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double theta_0, double phi_0)
{

  int ln_num = hierarchy->getNumberOfLevels();
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  patch_work_mpi_rank = -1, patch_work_level = -1, local_id = -1;

  patch_work_id = hier::GlobalId();

  for (int ln = 0; ln < ln_num; ln++)
  {
    boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    
    for (hier::PatchLevel::iterator p(level->begin());
         p != level->end(); ++p)
    {
      const boost::shared_ptr<hier::Patch>& patch = *p;
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      initPData(patch);
      initMDA(patch);
      
      const hier::Box& box = patch->getBox();
  
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];

      if(patch_work_i >= 0 &&
         (patch_work_i < lower[0] || patch_work_i > upper[0]
         || patch_work_j < lower[1] || patch_work_j > upper[1]
          || patch_work_k < lower[2] || patch_work_k > upper[2]))
      {
        continue;
      }
      
      const real_t * dx = &(patch_geometry->getDx())[0];

      boost::shared_ptr<pdat::CellData<real_t>> w_pdata (    
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  
          patch->getPatchData(w_idx)));
      arr_t w = pdat::ArrayDataAccess::access<DIM, double>(  
        w_pdata->getArrayData());

      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {

            if(w(i, j, k) > 0 && belowTheSurface(i, j, k))
            {
              real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0 - origin[0];
              real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0 - origin[1];
              real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0 - origin[2];

              real_t r = sqrt(pw2(x) + pw2(y) + pw2(z)) - F_a(i, j, k);

              real_t phi = atan(y / x);
              real_t theta = acos(z / sqrt(pw2(x) + pw2(y) + pw2(z))); 


              if( r * acos(cos(theta) * cos(theta_0) + sin(theta) * sin(theta_0) * cos(phi-phi_0)) < min_d + EPS)
              //if(fabs(sqrt(PW2(r *(theta - theta_0 ) ) + PW2(r*(phi-phi_0))) ) < min_d + EPS)
              //if(fabs(sqrt(PW2(r *(theta - theta_0 ) ) + PW2(r*(phi-phi_0))) - min_d) < EPS)
              {
                min_d = r * acos(cos(theta) * cos(theta_0) + sin(theta) * sin(theta_0) * cos(phi-phi_0));
                patch_work_i = i;
                patch_work_j = j;
                patch_work_k = k;
                patch_work_mpi_rank = mpi.getRank();
                patch_work_level = ln;
                local_id = patch->getLocalId().getValue();
              }
            }
          }
        }
      }

    }
  }

  
  mpi.Barrier();

  
  if (mpi.getSize() > 1)
  {
    mpi.AllReduce(&patch_work_mpi_rank, 1, MPI_MAX);
  }


  mpi.Barrier();

  if (mpi.getSize() > 1)
  {
    mpi.Bcast(&patch_work_i, 1, MPI_INT, patch_work_mpi_rank);
    mpi.Bcast(&patch_work_j, 1, MPI_INT, patch_work_mpi_rank);
    mpi.Bcast(&patch_work_k, 1, MPI_INT, patch_work_mpi_rank);
    mpi.Bcast(&patch_work_level, 1, MPI_INT, patch_work_mpi_rank);
    mpi.Bcast(&local_id, 1, MPI_INT, patch_work_mpi_rank);
  }
  mpi.Barrier();


  

}

void Horizon::initGridding(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  min_d = INF;
  real_t max_r = findMaxHorizonRadius(hierarchy, PI / 2.0, 0);

  if(is_sphere)
    tbox::pout<<"For spherical horizon, using the radius "<<max_r<<"\n";
  findPatch(hierarchy, PI / 2.0 , 0);
  
  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  int ln_num = hierarchy->getNumberOfLevels();

  real_t dx = tbox::MathUtilities<double>::Min(
    tbox::MathUtilities<double>::Min(grid_geometry.getDx()[0], grid_geometry.getDx()[1]),
    grid_geometry.getDx()[2]) / (real_t)(1<<(ln_num-1));
  
  n_theta = PI * max_r / dx;

  n_phi = 2.0 * PI * max_r / dx;

  tbox::pout<<"Dibviding the space into n_theta = "<<n_theta
            <<" and n_phi = "<<n_phi<<"\n";

  // initializing spherical mesh
  k_theta.resize(n_theta);
  k_phi.resize(n_theta);
  k_L.resize(n_theta);
  
  G111.resize(n_theta*2);
  G112.resize(n_theta*2);
  G122.resize(n_theta*2);
  G211.resize(n_theta*2);
  G212.resize(n_theta*2);
  G222.resize(n_theta*2);
  ah_radius.resize(n_theta*2);
  
  for(int i = 0; i < 2*n_theta; i++)
  {
    if(i < n_theta)
    {
      k_theta[i].resize(n_phi);
      k_phi[i].resize(n_phi);
      k_L[i].resize(n_phi);
    }
    G111[i].resize(n_phi*2);
    G112[i].resize(n_phi*2);
    G122[i].resize(n_phi*2);
    G211[i].resize(n_phi*2);
    G212[i].resize(n_phi*2);
    G222[i].resize(n_phi*2);
    ah_radius[i].resize(n_phi*2);
  }

}

bool solve_linear_eqn(int n,double a[][3],double b[])
{
  int i,j,k;
  double maxp;
  for (k = 0; k < n; k++)
  {
    maxp = a[k][k];

    for(i = k+1; i < n; i++)
    {
      double r = a[i][k] / maxp;
      for(j = k; j < n; j++)
        a[i][j] -= a[k][j] * r;
    }
  }

  
  for(i = n-1; i >= 0; i--)
  {
    double temp = 0;
    for(j = i+1; j < n; j++)
      temp += a[i][j] * b[j];
    if(fabs(a[i][i]) > EPS *10.0)
      b[i] = temp / a[i][i];
    else
      b[i] = 1.0;
  }
  
  return 1;
}
  
void Horizon::findM(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, double x[], BSSN *bssn)
{
  //  using namespace std::literals::complex_literals;

  transportKillingPhi(hierarchy, n_theta/2, n_phi, 1, 0, 0, bssn);

  Eigen::Matrix3d M;
  
  M(0, 0) = k_theta[n_theta/2][0];
  M(1, 0) = k_phi[n_theta/2][0];
  M(2, 0) = k_L[n_theta/2][0];


  transportKillingPhi(hierarchy, n_theta/2, n_phi, 0, 1, 0, bssn);

  M(0, 1) = k_theta[n_theta/2][0];
  M(1, 1) = k_phi[n_theta/2][0];
  M(2, 1) = k_L[n_theta/2][0];


  
  transportKillingPhi(hierarchy, n_theta/2, n_phi, 0, 0, 1, bssn);

  M(0, 2) = k_theta[n_theta/2][0]; 
  M(1, 2) = k_phi[n_theta/2][0];
  M(2, 2) = k_L[n_theta/2][0];



  tbox::pout<<"\n";
  tbox::pout << "Here is the matrix m:\n" << M << "\n";

  // solve the eigenvalue equation to get 3 eigenvalues;

  Eigen::EigenSolver<Eigen::Matrix3d> eigensolver(M);

  if (eigensolver.info() != Eigen::Success)
    TBOX_ERROR("Matrix does not have solution with identity eigenvalue!\n");

  Eigen::EigenSolver< Eigen::Matrix3d >::EigenvalueType e_val = eigensolver.eigenvalues();
  Eigen::EigenSolver< Eigen::Matrix3d >::EigenvectorsType e_vec = eigensolver.eigenvectors();

  tbox::pout<<"\n Eigenvalues are "<<e_val<<"\n";
  
  double dis_to_I = INF;
  int identity_idx=-1;
  for(int i = 0; i < 3; i++)
    if(pw2((e_val(i).real() - 1.0)) + pw2(e_val(i).imag()) < dis_to_I)
    {
      dis_to_I = pw2(e_val(i).real() - 1.0) + pw2(e_val(i).imag());
      identity_idx = i;
    }

  x[0] = e_vec(0, identity_idx).real();
  x[1] = e_vec(1, identity_idx).real();
  x[2] = e_vec(2, identity_idx).real();
  

  tbox::pout<<"Solution with identity eigenvalue is ("
            << x[0]<<" "<<x[1]<<" "<<x[2]<<")\n";
  
  
}

void Horizon::set_norm_values(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double theta, double phi, int theta_i, int phi_i, double r, KillingData *kd, BSSN * bssn)
{
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  cur_mpi_rank = -1;
  //if(patch_work_mpi_rank == mpi.getRank())
  int ln_num = hierarchy->getNumberOfLevels(), cur_mpi_level = -1, ln;

  real_t x = r * cos(phi) * sin(theta);
  real_t y = r * sin(phi) * sin(theta);
  real_t z = r * cos(theta);
        


  for (ln = ln_num - 1; ln >= 0; ln--)
  {
    boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    hier::PatchLevel::iterator p(level->begin());
    for (;p != level->end(); ++p)
    {
      const boost::shared_ptr<hier::Patch>& patch = *p;
      const hier::Box& box = patch->getBox();

      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];
      
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));


      const real_t * dx = &(patch_geometry->getDx())[0];
      int i0 = floor((x + origin[0] - domain_lower[0] ) / dx[0] - 0.5);
      int j0 = floor((y + origin[1] - domain_lower[1] ) / dx[1] - 0.5);
      int k0 = floor((z + origin[2] - domain_lower[2] ) / dx[2] - 0.5);

      if( i0 >= lower[0] && i0 <= upper[0]
          && j0 >= lower[1] && j0 <= upper[1]
          && k0 >= lower[2] && k0 <= upper[2])
      {
        cur_mpi_rank = mpi.getRank();

        cur_mpi_level = ln;
        
        initPData(patch);
        initMDA(patch);
        bssn->initPData(patch);
        bssn->initMDA(patch);

        BSSNData bd = {0};
        HorizonData hd = {0};

        for(int k = k0; k <= k0 + 1; k++)
        {
          real_t z0 = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0 - origin[2];

          COSMO_APPLY_TO_IJ_PERMS(HORIZON_DEFINE_TEMP_MJ);
          COSMO_APPLY_TO_IJ_PERMS(HORIZON_DEFINE_TEMP_KJ);

          HORIZON_DEFINE_TEMP_DFJ(1);
          HORIZON_DEFINE_TEMP_DFJ(2);
          HORIZON_DEFINE_TEMP_DFJ(3);
          
          for(int j = j0; j <= j0 + 1; j++)
          {
            real_t y0 = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0 - origin[1];

            COSMO_APPLY_TO_IJ_PERMS(HORIZON_DEFINE_TEMP_MI);
            COSMO_APPLY_TO_IJ_PERMS(HORIZON_DEFINE_TEMP_KI);
            HORIZON_DEFINE_TEMP_DFI(1);
            HORIZON_DEFINE_TEMP_DFI(2);
            HORIZON_DEFINE_TEMP_DFI(3);

            for(int i = i0; i <= i0 + 1; i++)
            {
              real_t x0 = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0 - origin[0];
              bssn->set_bd_values(i, j, k, &bd, dx);
              hd = getHorizonData(i, j, k, &bd, dx);
              COSMO_APPLY_TO_IJ_PERMS(HORIZON_INTERPOLATE_M_1);
              COSMO_APPLY_TO_IJ_PERMS(HORIZON_INTERPOLATE_K_1);
              HORIZON_INTERPOLATE_DF_1(1);
              HORIZON_INTERPOLATE_DF_1(2);
              HORIZON_INTERPOLATE_DF_1(3);
            }
      
            COSMO_APPLY_TO_IJ_PERMS(HORIZON_INTERPOLATE_M_2);
            COSMO_APPLY_TO_IJ_PERMS(HORIZON_INTERPOLATE_K_2);
            HORIZON_INTERPOLATE_DF_2(1);
            HORIZON_INTERPOLATE_DF_2(2);
            HORIZON_INTERPOLATE_DF_2(3);

          }
    
          COSMO_APPLY_TO_IJ_PERMS(HORIZON_INTERPOLATE_M_3);
          COSMO_APPLY_TO_IJ_PERMS(HORIZON_INTERPOLATE_K_3);
          HORIZON_INTERPOLATE_DF_3(1);
          HORIZON_INTERPOLATE_DF_3(2);
          HORIZON_INTERPOLATE_DF_3(3);
        }

        break;
      }
    }
    if(p != level->end())
      break;

  }

  real_t det = kd->m11 * kd->m22 * kd->m33 + kd->m12 * kd->m23 * kd->m13
    + kd->m12 * kd->m23 * kd->m13 - kd->m13 * kd->m22 * kd->m13
    - kd->m12 * kd->m12 * kd->m33 - kd->m23 * kd->m23 * kd->m11;
  
  kd->mi11 = (kd->m22 * kd->m33 - pw2(kd->m23)) / det;
  kd->mi22 = (kd->m11 * kd->m33 - pw2(kd->m13)) / det;
  kd->mi33 = (kd->m11 * kd->m22 - pw2(kd->m12)) / det;
  kd->mi12 = (kd->m13*kd->m23 - kd->m12*(kd->m33)) / det;
  kd->mi13 = (kd->m12*kd->m23 - kd->m13*(kd->m22)) / det;
  kd->mi23 = (kd->m12*kd->m13 - kd->m23*(kd->m11)) / det;

  
  mpi.Barrier();

  if (mpi.getSize() > 1)
  {
    mpi.AllReduce(&cur_mpi_level, 1, MPI_MAX);
  }
  mpi.Barrier();

  if(mpi.getSize() > 1 && ln != cur_mpi_level)
  {
    cur_mpi_rank = -1;
  }

  
  mpi.Barrier();

  if (mpi.getSize() > 1)
  {
    mpi.AllReduce(&cur_mpi_rank, 1, MPI_MAX);
  }
  mpi.Barrier();
  
  if(cur_mpi_rank == -1)
    TBOX_ERROR("Cannot find proper patch in set bd values\n");
  
}

real_t Horizon::interp_k_theta(
  double theta, double phi)
{
  if(theta > PI) theta = 2* PI - theta;
  if(theta < 0) theta = -theta;

  if(phi > 2.0 * PI) phi -= 2.0 * PI;
  if(phi < 0) phi = 2.0 * PI + phi;

  double d_theta = PI / (double)n_theta;
  double d_phi = 2.0 * PI / (double)n_phi;
  int theta_i0 = floor((theta) / d_theta - 0.5);
  int phi_i0 = floor((phi) / d_phi - 0.5);

  phi_i0 = (phi_i0 + n_phi)%n_phi;

  if(phi_i0 >= n_phi || phi_i0 < 0 || theta_i0 >= n_theta || theta_i0 < 0)
  {
    TBOX_ERROR("Step is too large!\n "<<theta<<" "<<phi<<" "<<theta_i0<<" "<<phi_i0<<"\n");
  }

  double res = 0;
  for(int theta_i = theta_i0; theta_i <= (theta_i0+1)%n_theta; theta_i = (theta_i +1 )%n_theta)
  {
    double theta0 = PI * ((double) theta_i +0.5) / (double)n_theta;

    double temp = 0;
    for(int phi_i = phi_i0; phi_i != (phi_i0+2)%n_phi; phi_i = (phi_i +1 )%n_phi)
    {
      double phi0 = 2.0 * PI * ((double) phi_i +0.5) / (double)n_phi;

      double phi_dist = phi0 - phi;
      if(phi_dist > d_phi + EPS)
        phi_dist -= 2.0 * PI;
      else if(phi_dist < -d_phi - EPS)
        phi_dist += 2.0 * PI;
      if(phi_dist > d_phi + EPS || phi_dist < -d_phi - EPS)
        TBOX_ERROR("Error in calculating phi_dist "<<phi0<<" "<<phi<<" "<<phi_dist<<" "<<d_phi<<"\n");
      temp += k_theta[theta_i][phi_i] *        
        ( (d_phi - (phi_dist ) * (2.0 * ( (phi_i - phi_i0 + n_phi)%n_phi ) -1.0)) / d_phi);
    }
    res += temp *
      ( (d_theta - (theta0 - theta ) * (2.0 * (theta_i - theta_i0) -1.0)) / d_theta);
  }
  
  return res;
}

real_t Horizon::interp_k_phi(
  double theta, double phi)
{
  if(theta > PI) theta = 2* PI - theta;
  if(theta < 0) theta = -theta;

  if(phi > 2.0 * PI) phi -= 2.0 * PI;
  if(phi < 0) phi = 2.0 * PI + phi;

  double d_theta = PI / (double)n_theta;
  double d_phi = 2.0 * PI / (double)n_phi;
  int theta_i0 = floor((theta) / d_theta - 0.5);
  int phi_i0 = floor((phi) / d_phi - 0.5);

  phi_i0 = (phi_i0 + n_phi)%n_phi;
  
  if(phi_i0 > n_phi || phi_i0 < 0 || theta_i0 > n_theta || theta_i0 < 0)
  {
    TBOX_ERROR("EEEE "<<theta<<" "<<phi<<" "<<theta_i0<<" "<<phi_i0<<"\n");
  }
  double res = 0;

  for(int theta_i = theta_i0; theta_i <= (theta_i0+1)%n_theta;  theta_i = (theta_i +1 )%n_theta)
  {
    double theta0 = PI * ((double) theta_i +0.5) / (double)n_theta;

    double temp = 0;
    for(int phi_i = phi_i0; phi_i != (phi_i0+2)%n_phi; phi_i = (phi_i +1 )%n_phi)
    {
      double phi0 = 2.0 * PI * ((double) phi_i +0.5) / (double)n_phi;

      double phi_dist = phi0 - phi;
      if(phi_dist > d_phi + EPS)
        phi_dist -= 2.0 * PI;
      else if(phi_dist < -d_phi - EPS)
        phi_dist += 2.0 * PI;
      if(phi_dist > d_phi + EPS || phi_dist < -d_phi - EPS)
        TBOX_ERROR("Error in calculating phi_dist "<<phi0<<" "<<phi<<" "<<phi_dist<<" "<<d_phi<<"\n");
      temp += k_phi[theta_i][phi_i] *        
        ( (d_phi - (phi_dist ) * (2.0 * ( (phi_i - phi_i0 + n_phi)%n_phi ) -1.0)) / d_phi);

    }
    res += temp *
      ( (d_theta - (theta0 - theta ) * (2.0 * (theta_i - theta_i0) -1.0)) / d_theta);
  }
  
  return res;

}

void Horizon::normKilling()
{
  double c = getNormFactor();

  tbox::pout<<"Norm factor "<<c<<"\n";
  for(int i = 0 ; i < n_theta; i++)
    for(int j = 0; j < n_phi; j ++)
    {
      k_theta[i][j] *= c;
      k_phi[i][j] *= c;
    }
}

real_t Horizon::getNormFactor()
{
  double theta = PI/2, phi = 0;

  real_t dtheta = PI / (double)n_theta;
  real_t pre_phi = 0;

  real_t max_abs = 0;

  for(int i = 0; i < n_theta; i ++)
    for(int j = 0; j < n_phi; j ++)
      max_abs = std::max(max_abs, std::max(fabs(k_phi[i][j]), fabs(k_theta[i][j])));
  
  real_t dt = 0.01 * dtheta / max_abs;

  tbox::pout<<"Time interval for process of getting norm factor is "
            <<dt<<"\n";

  
  real_t t = 0;
  // advance theta and phi
  while( (SIGN(pre_phi) * SIGN(phi) >= 0)
         && (SIGN(2.0 * PI - pre_phi) * SIGN(2.0 * PI - phi) >= 0 )
         &&(SIGN(-2.0 * PI - pre_phi) * SIGN(-2.0 * PI - phi) >= 0 )) 
  {
    pre_phi = phi;

    double theta_0 = theta;
    double phi_0 = phi;
    
    real_t k1_theta = interp_k_theta(theta, phi);
    real_t k1_phi = interp_k_phi(theta, phi);

    theta = theta_0 + k1_theta * dt / 2.0, phi = phi_0 + k1_phi * dt / 2.0;
    real_t k2_theta = interp_k_theta(theta , phi);
    real_t k2_phi = interp_k_phi(theta , phi);

    theta = theta_0 + k2_theta * dt / 2.0, phi = phi_0 + k2_phi * dt / 2.0;
    real_t k3_theta = interp_k_theta(theta, phi);
    real_t k3_phi = interp_k_phi(theta, phi);

    theta = theta_0 + k3_theta * dt, phi = phi_0 + k3_phi * dt;
    real_t k4_theta = interp_k_theta(theta , phi);
    real_t k4_phi = interp_k_phi(theta , phi);

    theta = theta_0 + dt * (k1_theta + 2.0 * k2_theta + 2.0 * k3_theta + k4_theta) / 6.0;
    phi = phi_0 + dt * (k1_phi + 2.0 * k2_phi + 2.0 * k3_phi + k4_phi) / 6.0;
    t += dt;
  }

  tbox::pout<<"When \phi gets back, the corresponding theta is "<<theta<<"\n";
  return t/ (2.0 * PI) ;
}

// here must be killing "vector" (not one form)
real_t Horizon::angularMomentum(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn)
{
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  real_t dtheta = PI / (double)n_theta;
  real_t dphi = 2.0 * PI / (double)n_phi;
  real_t res = 0;
  for(int theta_i = 0; theta_i < n_theta; theta_i++)
  {
    for(int phi_i = 0; phi_i < n_phi; phi_i++)
    {
      double theta = PI * ((double) theta_i +0.5) / (double)n_theta;
      double phi = 2.0 * PI * ((double) phi_i +0.5) / (double)n_phi;

      real_t st = sin(theta);
      real_t ct = cos(theta);
      real_t sp = sin(phi);
      real_t cp = cos(phi);
      
      real_t r = getRadius(theta_i*2, phi_i*2);
      
      KillingData kd = {0};

      // set_norm_values(hierarchy, theta, phi, theta_i, phi_i,
      //                 ah_radius[theta_i*2][phi_i*2], &kd, bssn);
      set_norm_values(hierarchy, theta, phi, theta_i, phi_i,
                      getRadius(theta_i*2, phi_i*2), &kd, bssn);
      

      real_t k1 = r * (k_theta[theta_i][phi_i] * ct * cp
                                                   - k_phi[theta_i][phi_i] * sp * st);
      real_t k2 = r * (k_theta[theta_i][phi_i] * ct * sp
                                                   + k_phi[theta_i][phi_i] * cp * st);
      real_t k3 = - r * k_theta[theta_i][phi_i] * st;

      
      real_t s1 = (kd.mi11 * kd.d1F + kd.mi12 * kd.d2F + kd.mi13 * kd.d3F)
        / (sqrt((kd.mi11 * kd.d1F * kd.d1F + kd.mi22 * kd.d2F * kd.d2F + kd.mi33 * kd.d3F *kd.d3F
                 + 2.0 * (kd.mi12 * kd.d1F * kd.d2F + kd.mi13 * kd.d1F * kd.d3F + kd.mi23 * kd.d2F * kd.d3F))));
      real_t s2 = (kd.mi21 * kd.d1F + kd.mi22 * kd.d2F + kd.mi23 * kd.d3F)
        / (sqrt((kd.mi11 * kd.d1F * kd.d1F + kd.mi22 * kd.d2F * kd.d2F + kd.mi33 * kd.d3F *kd.d3F
                 + 2.0 * (kd.mi12 * kd.d1F * kd.d2F + kd.mi13 * kd.d1F * kd.d3F + kd.mi23 * kd.d2F * kd.d3F))));
      real_t s3 = (kd.mi31 * kd.d1F + kd.mi32 * kd.d2F + kd.mi33 * kd.d3F)
        / (sqrt((kd.mi11 * kd.d1F * kd.d1F + kd.mi22 * kd.d2F * kd.d2F + kd.mi33 * kd.d3F *kd.d3F
                 + 2.0 * (kd.mi12 * kd.d1F * kd.d2F + kd.mi13 * kd.d1F * kd.d3F + kd.mi23 * kd.d2F * kd.d3F))));

      real_t K11 = kd.K11, K12 = kd.K12, K13 = kd.K13;
      real_t K22 = kd.K22, K23 = kd.K23, K33 = kd.K33;

      kd = {0};
      // set_kd_values(hierarchy, theta, phi, theta_i*2, phi_i*2,
      //               ah_radius[theta_i*2][phi_i*2], &kd, bssn);
      set_kd_values(hierarchy, theta, phi, theta_i*2, phi_i*2,
                    r, &kd, bssn);

      double det = kd.q11 * kd.q22 - kd.q12 * kd.q12;
      
      res += (k1 * s1 * K11 + k2 * s2 * K22 + k3 * s3 * K33
        + k1 * s2 * K12 + k1 * s3 * K13 + k2 * s3 * K23
              + k2 * s1 * K12 + k3 * s1 * K13 + k3 * s2 * K23) * sqrt(det) * dtheta * dphi;

      mpi.Barrier();
      if (mpi.getSize() > 1 ) {
        mpi.Bcast(&res, 1, MPI_DOUBLE, cur_mpi_rank);
      }
      mpi.Barrier();
    }
  }
  return res / 8.0 / PI;
}

void Horizon::convertToVector(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn)
{
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  for(int theta_i = 0; theta_i < n_theta; theta_i++)
  {
    for(int phi_i = 0; phi_i < n_phi; phi_i++)
    {
      double k_theta0 = k_theta[theta_i][phi_i];
      double k_phi0 = k_phi[theta_i][phi_i];
      double theta = PI * ((double) theta_i +0.5) / (double)n_theta;
      double phi = 2.0 * PI * ((double) phi_i +0.5) / (double)n_phi;

      KillingData kd = {0};
      //      set_kd_values(hierarchy, theta, phi, theta_i*2, phi_i*2, ah_radius[theta_i*2][phi_i*2], &kd, bssn);
      set_kd_values(hierarchy, theta, phi, theta_i*2, phi_i*2, getRadius(theta_i*2, phi_i*2), &kd, bssn);
      
      k_theta[theta_i][phi_i] = kd.qi11 * k_theta0 + kd.qi12 * k_phi0;
      k_phi[theta_i][phi_i] = kd.qi12 * k_theta0 + kd.qi22 * k_phi0;

      mpi.Barrier();
      if (mpi.getSize() > 1 ) {
        mpi.Bcast(&k_theta[theta_i][phi_i], 1, MPI_DOUBLE, cur_mpi_rank);
        mpi.Bcast(&k_phi[theta_i][phi_i], 1, MPI_DOUBLE, cur_mpi_rank);
      }
      mpi.Barrier();

    }
  }


}

real_t Horizon::area(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn)
{
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  real_t dtheta = PI / (double)n_theta;
  real_t dphi = 2.0 * PI / (double)n_phi;
  real_t res = 0;
  for(int theta_i = 0; theta_i < n_theta; theta_i++)
  {
    for(int phi_i = 0; phi_i < n_phi; phi_i++)
    {
      double theta = PI * ((double) theta_i +0.5) / (double)n_theta;
      double phi = 2.0 * PI * ((double) phi_i +0.5) / (double)n_phi;

      real_t st = sin(theta);
      real_t ct = cos(theta);
      real_t sp = sin(phi);
      real_t cp = cos(phi);
      
      real_t r = getRadius(theta_i*2, phi_i*2);
      
      KillingData kd = {0};

      set_kd_values(hierarchy, theta, phi, theta_i*2, phi_i*2,
                    r, &kd, bssn);

      double det = kd.q11 * kd.q22 - kd.q12 * kd.q12;
      
      res += sqrt(det) * dtheta * dphi;

      mpi.Barrier();
      if (mpi.getSize() > 1 ) {
        mpi.Bcast(&res, 1, MPI_DOUBLE, cur_mpi_rank);
      }
      mpi.Barrier();
    }
  }
  return res;

}


void Horizon::findKilling(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn)
{
  tbox::pout<<"Starting the process of finding Killing vectors!\n";

  initGridding(hierarchy);

  initG(hierarchy, bssn);

  double x[3] = {0};
  // Finding transport matrix and initializing eigen vector
  findM(hierarchy, x, bssn);
  
  transportKillingPhi(hierarchy, n_theta/2, n_phi - 1, x[0], x[1], x[2], bssn);

  for(int i = 0; i < n_phi; i++)
  {
    transportKillingTheta(
      hierarchy, i, k_theta[n_theta/2][i], k_phi[n_theta/2][i], k_L[n_theta/2][i], bssn);
  }

  
  convertToVector(hierarchy, bssn);

  normKilling();

  double angular_m = angularMomentum(hierarchy, bssn);
  tbox::pout<<"Angular momentum is "<<angular_m<<"\n";

  
  double a = area(hierarchy, bssn);

  double R_Delta = sqrt(a / (4.0 * PI));
  
  double mass = sqrt(pw2(pw2(R_Delta)) + 4.0 * pw2(angular_m)) / (2.0 * R_Delta);

  tbox::pout<<"Mass is "<<mass<<"\n";
}
  
}
