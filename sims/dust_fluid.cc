#include "dust_fluid.h"
#include "../cosmo_includes.h"
#include "../utils/math.h"
#include "SAMRAI/xfer/PatchLevelBorderFillPattern.h"

using namespace SAMRAI;

namespace cosmo
{
/**
 * @brief Constructing CosmoSim object
 * 
 * @param hierarchy 
 * @param dimenstion
 * @param input database, which includes ALL sub databases
 * @param IO stream
 * @param name of simulation type
 * @param name of output file for VisIt
 */

DustFluidSim::DustFluidSim(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  const tbox::Dimension& dim_in,
  std::shared_ptr<tbox::InputDatabase>& input_db_in,
  std::ostream* l_stream_in,
  std::string simulation_type_in,
  std::string vis_filename_in):CosmoSim(
    hierarchy,
    dim_in, input_db_in, l_stream_in, simulation_type_in, vis_filename_in),
    cosmo_dust_fluid_db(input_db_in->getDatabase("DustFluidSim")),
    freeze_fluid_step(cosmo_dust_fluid_db->getIntegerWithDefault("freeze_fluid_step", 999999999))
{
  if(!USE_DUST_FLUID)
    TBOX_ERROR("Macro USE_DUST_FLUID needs to be true in order to evolve dust fluid!");
  t_init->start();

  std::string bd_type = cosmo_dust_fluid_db->getString("boundary_type");
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  
  if(bd_type == "periodic")
  {
    cosmoPS = new periodicBD(dim, bd_type);
  }
  else
    TBOX_ERROR("Unsupported boundary type!\n");


  cosmo_dust_fluid_db = input_db->getDatabase("DustFluidSim");


  tbox::pout<<"Running 'dust' type simulation.\n";

    gradient_indicator_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable(gradient_indicator), variable_db->getContext("ACTIVE"));

  
  // adding all fields to a list
  bssnSim->addFieldsToList(variable_id_list);

  dustFluidSim = new DustFluid(
    hierarchy, dim, input_db->getDatabase("DustFluid"), lstream);

  dustFluidSim->addFieldsToList(variable_id_list);
  //  variable_id_list.push_back(dustFluidSim->DIFFD_a_idx);
  
  variable_id_list.push_back(weight_idx);

  hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

  tbox::RestartManager::getManager()->registerRestartItem(simulation_type_in,
                                                          this);

  for(int i = 0; i < static_cast<idx_t>(variable_id_list.size()); i++)
  {
    hier::PatchDataRestartManager::getManager()->
      registerPatchDataForRestart(variable_id_list[i]);
  }

  if(tbox::RestartManager::getManager()->isFromRestart())
    getFromRestart();

  t_init->stop();  
}

DustFluidSim::~DustFluidSim() {
}

  
void DustFluidSim::init()
{
  
}

/**
 * @brief      Set dust initial conditions
 *
 */
void DustFluidSim::setICs(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  tbox::plog<<"Setting initial conditions (ICs).";

  TBOX_ASSERT(gridding_algorithm);
  
  gridding_algorithm->printClassData(tbox::plog);

  bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();

  if(is_from_restart)
    hierarchy->initializeHierarchy();
  
  // set initial condition by calling function initializeLevelData()
  gridding_algorithm->makeCoarsestLevel(cur_t);

  if(is_from_restart)
  {
    std::vector<int> tag_buffer(hierarchy->getMaxNumberOfLevels());
    for (idx_t ln = 0; ln < static_cast<int>(tag_buffer.size()); ++ln) {
      tag_buffer[ln] = 1;
    }
    gridding_algorithm->regridAllFinerLevels(
      0,
      tag_buffer,
      0,
      cur_t);
  }

  
  // regrid initial hierarchy if needed
  while(!is_from_restart &&
        hierarchy->getNumberOfLevels() < hierarchy->getMaxNumberOfLevels())
  {
    int pre_level_num = hierarchy->getNumberOfLevels();
    std::vector<int> tag_buffer(hierarchy->getMaxNumberOfLevels());
    for (idx_t ln = 0; ln < static_cast<int>(tag_buffer.size()); ++ln) {
      tag_buffer[ln] = 1;
    }
    gridding_algorithm->regridAllFinerLevels(
      0,
      tag_buffer,
      0,
      cur_t);

    int post_level_num = hierarchy->getNumberOfLevels();
    // no new level is created
    if(post_level_num == pre_level_num) break;
  }
#if USE_COSMOTRACE
  //  ray->initAll(hierarchy, bssnSim, dustFluidSim);
#endif

  /* Set vector weight. */
  computeVectorWeights(hierarchy);
  
  tbox::plog<<"Finished setting ICs. with hierarchy has "
            <<hierarchy->getNumberOfLevels()<<" levels\n";
}

void DustFluidSim::initDustFluidStep(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  bssnSim->stepInit(hierarchy);
  bssnSim->clearSrc(hierarchy);
  dustFluidSim->addBSSNSrc(bssnSim, hierarchy);
  dustFluidSim->addDerivedFields(bssnSim, hierarchy);

#if USE_COSMOTRACE
  ray->clearParticlesLivingInGhostCells(hierarchy);
  if(ray_insert_step.empty() && step == 0)
    ray->initAll(hierarchy, bssnSim, dustFluidSim);
  else
  {
    for(int i = 0; i < ray_insert_step.size(); i ++)
      if(step == ray_insert_step[i])
      {
        ray->initAll(hierarchy, bssnSim, dustFluidSim);
        break;
      }        
  }
#endif
}

/**
 * @brief      initilize newly created level
 *             set value directly if it's possible and return true
 *             do nothing when it's not possible and return false
 * @param hierarchy
 * @param level index
 */
bool DustFluidSim::initLevel(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln)
{
  std::string ic_type = cosmo_dust_fluid_db->getString("ic_type");

  std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

  math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);

  // zero all fields
  bssnSim->clearField(hierarchy, ln);
  bssnSim->clearSrc(hierarchy, ln);
  bssnSim->clearGen1(hierarchy, ln);

  hcellmath.setToScalar(dustFluidSim->DF_D_a_idx, 0, 0);
  hcellmath.setToScalar(dustFluidSim->DF_S1_a_idx, 0, 0);
  hcellmath.setToScalar(dustFluidSim->DF_S2_a_idx, 0, 0);
  hcellmath.setToScalar(dustFluidSim->DF_S3_a_idx, 0, 0);
  hcellmath.setToScalar(dustFluidSim->DF_E_a_idx, 0, 0);
  
  if(ic_type == "fluid_for_BHL")
  {
    if(ln > 0) return false;
    dustFluidSim->dust_fluid_ic_set_fluid_for_BHL(hierarchy,ln, input_db->getDatabase("DustFluidSim"));
    return true;
  }
  else
    TBOX_ERROR("Undefined IC type!\n");

  return false;
}


/**
 * @brief initializeLevelData when there is new level created
 * 
 * @param Hierarchy to initialize 
 * @param level index
 * @param the time to initialize the level
 * @param whether the level can be refined
 * @param whether level is being introduced for the first time
 * @param level to copy data from
 * @param whether the level has been alocated memories
 */
  
void DustFluidSim::initializeLevelData(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  const int ln,
  const double init_data_time,
  const bool can_be_refined,
  const bool initial_time,
  const std::shared_ptr<hier::PatchLevel>& old_level,
  const bool allocate_data)
{
   NULL_USE(can_be_refined);
   NULL_USE(initial_time);

   std::shared_ptr<hier::PatchHierarchy> patch_hierarchy(hierarchy);

   /*
    * Reference the level object with the given index from the hierarchy.
    */
   std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
   

   math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);
      
   
   if (allocate_data)
   {
     bssnSim->allocField(patch_hierarchy, ln);
     bssnSim->allocSrc(patch_hierarchy, ln);
     bssnSim->allocGen1(patch_hierarchy, ln);
     dustFluidSim->alloc(patch_hierarchy, ln);
     level->allocatePatchData(weight_idx);
#if USE_COSMOTRACE
     ray->allocParticles(patch_hierarchy, ln);
#endif

   }

   // marks whether we have solved initial value for certain level,
   // if so, there is no need to interpolate from coarser level
   bool has_initial = false;
   
   //at beginning, initialize new level
   if(fabs(init_data_time - starting_t)< EPS)
   {
     if(step != starting_step)
       TBOX_ERROR("Level is initialized after 0 step!");
     has_initial = initLevel(patch_hierarchy, ln);
   }
   bssnSim->clearSrc(patch_hierarchy, ln);
   bssnSim->clearGen1(patch_hierarchy, ln);
   dustFluidSim->clearDerivedFields(patch_hierarchy, ln);
   /*
    * Refine solution data from coarser level and, if provided, old level.
    */
   
   xfer::RefineAlgorithm refiner;

   std::shared_ptr<hier::RefineOperator> accurate_refine_op =
     space_refine_op;
     
   TBOX_ASSERT(accurate_refine_op);

   //registering refine variables
   //BSSN_APPLY_TO_FIELDS_ARGS(DUS_REGISTER_SPACE_REFINE_A,refiner,accurate_refine_op);

   bssnSim->registerRKRefinerActive(refiner, accurate_refine_op);

#if USE_COSMOTRACE
   ray->registerRKRefiner(refiner, particle_refine_op);
#endif

   dustFluidSim->registerRKRefinerActive(refiner, accurate_refine_op);

#if USE_COSMOTRACE
   ray->registerRKRefiner(refiner, particle_refine_op);
#endif

   
   std::shared_ptr<xfer::RefineSchedule> refine_schedule;

   level->getBoxLevel()->getMPI().Barrier();
   if (ln > 0 && (!has_initial))
   {
     /*
      * Include coarser levels in setting data
      */
     refine_schedule =
       refiner.createSchedule(
         level,
         old_level,
         ln - 1,
         hierarchy,
         cosmoPS);
   }
   else
   {
     /*
      * There is no coarser level, and source data comes only
      * from old_level, if any.
      */
     if (old_level)
     {
       refine_schedule =
         refiner.createSchedule(level,
                                old_level,
                                NULL);
     }
     else
     {
       refine_schedule =
         refiner.createSchedule(level,
                                level,
                                NULL);
     }
   }
   level->getBoxLevel()->getMPI().Barrier();
     

   if (refine_schedule)
   {
     refine_schedule->fillData(0.0);
     // It is null if this is the bottom level.
   }
   else
   {
     TBOX_ERROR(
       "Can not get refine schedule, check your code!\n");
   }
 
   bssnSim->copyAToP(hcellmath);
   dustFluidSim->copyAToP(hcellmath);   
   level->getBoxLevel()->getMPI().Barrier();
}

/**
 * @brief tag the grids that need to be refined
 *
 */  
void DustFluidSim::applyGradientDetector(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy_,
   const int ln,
   const double error_data_time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation)
{
  NULL_USE(uses_richardson_extrapolation);
  NULL_USE(error_data_time);
  NULL_USE(initial_time);

  if (lstream) {
    *lstream
      << "VaccumSim("  << ")::applyGradientDetector"
      << std::endl;
  }
  hier::PatchHierarchy& hierarchy = *hierarchy_;
  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy.getGridGeometry()));
  double max_der_norm = 0;
  hier::PatchLevel& level =
    (hier::PatchLevel &) * hierarchy.getPatchLevel(ln);
  int ntag = 0, ntotal = 0;
  //double maxestimate = 0;
  for (hier::PatchLevel::iterator pi(level.begin());
       pi != level.end(); ++pi)
  {
    const std::shared_ptr<hier::Patch> & patch = *pi;

    const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry()));


    std::shared_ptr<pdat::CellData<real_t> > f_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<real_t>, hier::PatchData>(
        patch->getPatchData(gradient_indicator_idx)));

    arr_t f =
      pdat::ArrayDataAccess::access<DIM, real_t>(
        f_pdata->getArrayData());
    std::shared_ptr<pdat::CellData<int> > tag_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
        patch->getPatchData(tag_index)));
      
    MDA_Access<int, DIM, MDA_OrderColMajor<DIM>>  tag =
      pdat::ArrayDataAccess::access<DIM, int>(
        tag_pdata->getArrayData());

    ntotal += patch->getBox().numberCells().getProduct();

    const hier::Box& box = patch->getBox();
    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

#pragma omp parallel for collapse(2) reduction(+:ntag) reduction( max: max_der_norm)
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          tag(i, j, k) = 0;
          max_der_norm = tbox::MathUtilities<double>::Max(
            max_der_norm,
            derivative_norm(i, j, k, f));

          if(derivative_norm(i, j, k, f) > adaption_threshold )
          {
            tag(i, j, k) = 1;
            ++ntag;
          }
            
        }
      }
    }

  }
  const tbox::SAMRAI_MPI& mpi(hierarchy.getMPI());
  if (mpi.getSize() > 1)
  {
    mpi.AllReduce(&max_der_norm, 1, MPI_MAX);
    mpi.AllReduce(&ntag, 1, MPI_SUM);
    mpi.AllReduce(&ntotal, 1, MPI_SUM);
  }

  tbox::plog << "Adaption threshold is " << adaption_threshold << "\n";
  tbox::plog << "Number of cells tagged on level " << ln << " is "
             << ntag << "/" << ntotal << "\n";
  tbox::plog << "Max norm is " << max_der_norm << "\n";

}
void DustFluidSim::outputDustFluidStep(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{

  std::shared_ptr<appu::VisItDataWriter> visit_writer(
    new appu::VisItDataWriter(
      dim, "VisIt Writer", vis_filename + ".visit"));

  tbox::pout<<"step: "<<step<<"/"<<num_steps<<"\n";

#if CAL_WEYL_SCALS
  if(calculate_Weyl_scalars)
    bssnSim->cal_Weyl_scalars(hierarchy, weight_idx);
#endif
#if USE_COSMOTRACE
  //    ray->printAll(hierarchy, ray->pc_idx);
#endif
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
  bssnSim->output_L2_H_constaint(
    hierarchy, weight_idx, cosmoPS, 0.5);
  //  bssnSim->output_max_H_constaint(hierarchy, weight_idx);
 
  cosmo_io->registerVariablesWithPlotter(*visit_writer, step);
#if !USE_COSMOTRACE  
  cosmo_io->dumpData(hierarchy, *visit_writer, step, cur_t);
#else
  cosmo_io->dumpData(hierarchy, *visit_writer, step, cur_t, ray, vis_filename);
#endif
  cosmo_statistic->output_expansion_info(
      hierarchy,
      bssnSim, weight_idx, step, cur_t, max_horizon_radius);

  cosmo_statistic->output_conformal_avg(
    hierarchy,
    bssnSim, weight_idx, step, cur_t, max_horizon_radius);
  dustFluidSim->printWConstraint(bssnSim, hierarchy, weight_idx);
}

/**
 * @brief advance hierarchy from time "from_t" to "to_t"
 * 
 * @param hierarchy
 * @param starting time
 * @param ending time
 */  
void DustFluidSim::runDustFluidStep(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double from_t, double to_t)
{
  t_RK_steps->start();
    // Full RK step minus init()
  advanceLevel(hierarchy,
               0,
               from_t,
               to_t);

  t_RK_steps->stop();
}

/**
 * @brief  get dt for each step, currently will return the same value
 *
 */  
double DustFluidSim::getDt(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  return tbox::MathUtilities<double>::Min(
    tbox::MathUtilities<double>::Min(grid_geometry.getDx()[0], grid_geometry.getDx()[1]),
    grid_geometry.getDx()[2])
    * dt_frac;
  return (grid_geometry.getDx()[0]) * dt_frac;
}


/**
 * @brief run each step
 *
 */  
void DustFluidSim::runStep(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  runCommonStepTasks(hierarchy);
  double dt = getDt(hierarchy);
  initDustFluidStep(hierarchy);
  
  outputDustFluidStep(hierarchy);
  
  runDustFluidStep(hierarchy, cur_t, cur_t + dt);
  cur_t += dt;
}


void DustFluidSim::addBSSNExtras(
  const std::shared_ptr<hier::PatchLevel> & level)
{
  return;
}

void DustFluidSim::addBSSNExtras(
  const std::shared_ptr<hier::Patch> & patch)
{
  return;
}


void DustFluidSim::RKEvolvePatch(
  const std::shared_ptr<hier::Patch> & patch, real_t dt)
{
  bssnSim->initPData(patch);
  bssnSim->initMDA(patch);
  dustFluidSim->initPData(patch);
  dustFluidSim->initMDA(patch);
  const hier::Box& box = patch->getBox();
  
  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
  const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(  
    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  //initialize dx for each patch
  const real_t * dx = &(patch_geom->getDx())[0];


  if(step <= freeze_fluid_step)
  {
#pragma omp parallel for collapse(2)
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          BSSNData bd = {0};
          DustFluidData dd = {0};
          bssnSim->RKEvolvePt(i, j, k, bd, dx, dt);
          dustFluidSim->RKEvolvePt(i, j, k, bd, dd, dx, dt);
        }
      }
    }
  }
  else
  {
#pragma omp parallel for collapse(2)
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          BSSNData bd = {0};
          DustFluidData dd = {0};
          bssnSim->RKEvolvePt(i, j, k, bd, dx, dt);
        }
      }
    }    
  }
}

void DustFluidSim::RKEvolvePatchBD(
  const std::shared_ptr<hier::Patch> & patch, real_t dt)
{
    std::shared_ptr<hier::PatchGeometry> geom (patch->getPatchGeometry());

  idx_t codim = 1;

  // getting all codimension 1 boxes
  const std::vector<hier::BoundaryBox> & codim1_boxes =
    geom->getCodimensionBoundaries(codim);

  const idx_t n_codim1_boxes = static_cast<idx_t>(codim1_boxes.size());

  // if it has no codimention 1 boundary, it has no other type of boundaries
  if(n_codim1_boxes == 0) return;

  bssnSim->initPData(patch);
  bssnSim->initMDA(patch);
  dustFluidSim->initPData(patch);
  dustFluidSim->initMDA(patch);

  const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom( 
    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

    
  const hier::Box& patch_box = patch->getBox();


  for(int l = 0 ; l < n_codim1_boxes; l++)
  {
    hier::Box boundary_fill_box =
      geom->getBoundaryFillBox(
        codim1_boxes[l], patch_box, bssnSim->DIFFchi_a_pdata->getGhostCellWidth());

    if(boundary_fill_box.empty()) continue;


    idx_t l_idx = codim1_boxes[l].getLocationIndex();
    

    boundary_fill_box.shift(
      (hier::Box::dir_t)l_idx/2,
      (l_idx%2)?(-STENCIL_ORDER_WIDTH):STENCIL_ORDER_WIDTH);


    boundary_fill_box *= patch_box;

    //initialize dx for each patch
    const real_t * dx = &(patch_geom->getDx())[0];

    const idx_t * lower = &boundary_fill_box.lower()[0];
    const idx_t * upper = &boundary_fill_box.upper()[0];

    #pragma omp parallel for collapse(2)
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          BSSNData bd = {0};
          DustFluidData dd = {0};
          bssnSim->RKEvolvePtBd(i, j, k, bd, dx, dt, l_idx, codim);
          dustFluidSim->RKEvolvePtBd(i, j, k, bd, dd, dx, dt, l_idx, codim);
        }
      }
    }
  }
  /************************updating codim = 2 boundaries****************/
  codim = 2;

  const std::vector<hier::BoundaryBox> & codim2_boxes =
    geom->getCodimensionBoundaries(codim);

  const idx_t n_codim2_boxes = static_cast<idx_t>(codim2_boxes.size());




  for(int l = 0 ; l < n_codim2_boxes; l++)
  {
    hier::Box  boundary_fill_box =
      geom->getBoundaryFillBox(
        codim2_boxes[l], patch_box, bssnSim->DIFFchi_a_pdata->getGhostCellWidth());
    
    if(boundary_fill_box.empty()) continue;  


    idx_t l_idx = codim2_boxes[l].getLocationIndex();
    

    std::vector<idx_t> shift_vec;
    if(l_idx == 0 || l_idx == 2 || l_idx == 4 || l_idx == 6)
      shift_vec.push_back(STENCIL_ORDER_WIDTH);
    else if(l_idx == 1 || l_idx == 3 || l_idx == 5 || l_idx == 7)
      shift_vec.push_back(-STENCIL_ORDER_WIDTH);
    else
      shift_vec.push_back(0);
    
    if(l_idx == 0 || l_idx == 1 || l_idx == 8 || l_idx == 10)
      shift_vec.push_back(STENCIL_ORDER_WIDTH);
    else if(l_idx == 2 || l_idx == 3 || l_idx ==9 || l_idx == 11)
     shift_vec.push_back(-STENCIL_ORDER_WIDTH);
    else
      shift_vec.push_back(0);

    if( l_idx == 4 || l_idx == 5 || l_idx == 8 || l_idx == 9)
      shift_vec.push_back(STENCIL_ORDER_WIDTH);
    else if(l_idx == 6 || l_idx == 7 || l_idx == 10 || l_idx == 11)
      shift_vec.push_back(-STENCIL_ORDER_WIDTH);
    else
      shift_vec.push_back(0);

    boundary_fill_box.shift(hier::IntVector(shift_vec));

    boundary_fill_box *= patch_box;

    //initialize dx for each patch
    const real_t * dx = &(patch_geom->getDx())[0];

    const idx_t * lower = &boundary_fill_box.lower()[0];
    const idx_t * upper = &boundary_fill_box.upper()[0];

#pragma omp parallel for collapse(2)    
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          BSSNData bd = {0};
          DustFluidData dd = {0};
          bssnSim->RKEvolvePtBd(i, j, k, bd, dx, dt, l_idx, codim);
          dustFluidSim->RKEvolvePtBd(i, j, k, bd, dd, dx, dt, l_idx, codim);
        }
      }
    }
  }

  /************************updating codim = 3 boundaries****************/
  codim = 3;

  const std::vector<hier::BoundaryBox> & codim3_boxes =
    geom->getCodimensionBoundaries(codim);

  const idx_t n_codim3_boxes = static_cast<idx_t>(codim3_boxes.size());




  for(int l = 0 ; l < n_codim3_boxes; l++)
  {
    hier::Box boundary_fill_box =
      geom->getBoundaryFillBox(
        codim3_boxes[l], patch_box, bssnSim->DIFFchi_a_pdata->getGhostCellWidth());

    if(boundary_fill_box.empty()) continue;

    idx_t l_idx = codim3_boxes[l].getLocationIndex();

    std::vector<idx_t> shift_vec;
    
    if(l_idx == 0 || l_idx == 2 || l_idx == 4 || l_idx == 6)
      shift_vec.push_back(STENCIL_ORDER_WIDTH);
    else
      shift_vec.push_back(-STENCIL_ORDER_WIDTH);
    
    if(l_idx == 0 || l_idx == 1 || l_idx == 4 || l_idx == 5)
      shift_vec.push_back(STENCIL_ORDER_WIDTH);
    else
      shift_vec.push_back(-STENCIL_ORDER_WIDTH);

    if( l_idx == 0 || l_idx == 1 || l_idx == 2 || l_idx == 3)
      shift_vec.push_back(STENCIL_ORDER_WIDTH);
    else
      shift_vec.push_back(-STENCIL_ORDER_WIDTH);

    boundary_fill_box.shift(hier::IntVector(shift_vec));

    boundary_fill_box *= patch_box;

    const idx_t * lower = &boundary_fill_box.lower()[0];
    const idx_t * upper = &boundary_fill_box.upper()[0];

    //initialize dx for each patch
    const real_t * dx = &(patch_geom->getDx())[0];
    #pragma omp parallel for collapse(2)
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          BSSNData bd = {0};
          DustFluidData dd = {0};
          bssnSim->RKEvolvePtBd(i, j, k, bd, dx, dt, l_idx, codim);
          dustFluidSim->RKEvolvePtBd(i, j, k, bd, dd, dx, dt, l_idx, codim);
                    
        }
      }
    }
  }
}

/**
 * @brief RK evolve level
 * 
 * @param hierarchy
 * @param level index
 * @param starting time
 * @param ending time
 */  
void DustFluidSim::RKEvolveLevel(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln,
  double from_t,
  double to_t)
{
  const std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  const std::shared_ptr<hier::PatchLevel> coarser_level(
    ((ln>0)?(hierarchy->getPatchLevel(ln-1)):NULL));
  
  
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      bssnSim->prepareForK1(coarser_level, to_t);
      dustFluidSim->prepareForK1(coarser_level, to_t);
    }
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    //Evolve inner grids
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      RKEvolvePatch(patch, to_t - from_t);
    }

    // Evolve physical boundary
    // would not do anything if boundary is time independent

#if USE_COSMOTRACE
    ray->RKEvolvePatch(patch, bssnSim, dustFluidSim,to_t - from_t);
#endif

#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      RKEvolvePatchBD(patch, to_t - from_t);
    }
  }

  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
  pre_refine_schedules[ln]->fillData(to_t);

  bssnSim->clearSrc(hierarchy, ln);
  dustFluidSim->clearDerivedFields(hierarchy, ln);

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    bssnSim->K1FinalizePatch(patch);    
    addBSSNExtras(patch);
    if(step <= freeze_fluid_step)
    {
      dustFluidSim->K1FinalizePatch(patch);
      dustFluidSim->addDerivedFields(bssnSim, patch);
      dustFluidSim->addBSSNSrc(bssnSim,patch, false);
    }
#if USE_COSMOTRACE
    ray->K1FinalizePatch(patch);
#endif
  }

  bssnSim->set_norm(level);
  
  /**************Starting K2 *********************************/
  #if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      bssnSim->prepareForK2(coarser_level, to_t);
      dustFluidSim->prepareForK2(coarser_level, to_t);
    }
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    //Evolve inner grids
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      RKEvolvePatch(patch, to_t - from_t);
    }
    //Evolve physical boundary
    // would not do anything if boundary is time dependent
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      RKEvolvePatchBD(patch, to_t - from_t);
    }

#if USE_COSMOTRACE
    ray->RKEvolvePatch(patch, bssnSim, dustFluidSim, to_t - from_t);
#endif

  }

  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
  #if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
  pre_refine_schedules[ln]->fillData(to_t);
  bssnSim->clearSrc(hierarchy, ln);
  dustFluidSim->clearDerivedFields(hierarchy, ln);

    
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      bssnSim->K2FinalizePatch(patch);
      addBSSNExtras(patch);
      if(step <= freeze_fluid_step)
      {
        dustFluidSim->K2FinalizePatch(patch);
        dustFluidSim->addDerivedFields(bssnSim, patch);
        dustFluidSim->addBSSNSrc(bssnSim,patch, false);
      }

    }
#if USE_COSMOTRACE
    ray->K2FinalizePatch(patch);
#endif
  }

  bssnSim->set_norm(level);

  /**************Starting K3 *********************************/
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      bssnSim->prepareForK3(coarser_level, to_t);
      dustFluidSim->prepareForK3(coarser_level, to_t);
    }
    

  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      RKEvolvePatch(patch, to_t - from_t);
    }
    //Evolve physical boundary
    // would not do anything if boundary is time dependent
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      RKEvolvePatchBD(patch, to_t - from_t);
    }
#if USE_COSMOTRACE
    ray->RKEvolvePatch(patch, bssnSim, dustFluidSim, to_t - from_t);
#endif

  }

  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
  #if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
  pre_refine_schedules[ln]->fillData(to_t);
  bssnSim->clearSrc(hierarchy, ln);
  dustFluidSim->clearDerivedFields(hierarchy, ln);


  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      bssnSim->K3FinalizePatch(patch);
      addBSSNExtras(patch);
      if(step <= freeze_fluid_step)
      {
        dustFluidSim->K3FinalizePatch(patch);
        dustFluidSim->addDerivedFields(bssnSim, patch);
        dustFluidSim->addBSSNSrc(bssnSim,patch, false);
      }
    }
#if USE_COSMOTRACE
    ray->K3FinalizePatch(patch);
#endif
  }

  bssnSim->set_norm(level);

  /**************Starting K4 *********************************/
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      bssnSim->prepareForK4(coarser_level, to_t);
      dustFluidSim->prepareForK4(coarser_level, to_t);
    }

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      RKEvolvePatch(patch, to_t - from_t);
    }
    //Evolve physical boundary
    // would not do anything if boundary is time dependent
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      RKEvolvePatchBD(patch, to_t - from_t);
    }
#if USE_COSMOTRACE
    ray->RKEvolvePatch(patch, bssnSim, dustFluidSim, to_t - from_t);
#endif

  }

  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
  pre_refine_schedules[ln]->fillData(to_t);
  bssnSim->clearSrc(hierarchy, ln);
  dustFluidSim->clearDerivedFields(hierarchy, ln);


  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
#if USE_COSMOTRACE
    if(freeze_time_evolution == false)
#endif
    {
      bssnSim->K4FinalizePatch(patch);
      addBSSNExtras(patch);
      if(step <= freeze_fluid_step)
      {
        dustFluidSim->K4FinalizePatch(patch);
        dustFluidSim->addDerivedFields(bssnSim, patch);
        dustFluidSim->addBSSNSrc(bssnSim,patch, false);
      }
    }
#if USE_COSMOTRACE
    ray->K4FinalizePatch(patch);
#endif
  }

  bssnSim->set_norm(level);

}

/**
 * @brief advance a single level by:
 *        1. RK evolve this level(including evolve the interior and interpolate the boundary)
 *        2. evolve its son levels recursively
 *        3. doing coarsen operation
 */ 
void DustFluidSim::advanceLevel(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  int ln,
  double from_t,
  double to_t)
{
  if( ln >= hierarchy->getNumberOfLevels())
    return;
  
  //double dt = to_t - from_t;
  const std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  //updating extra fields before advancing any level
  addBSSNExtras(level);

  bssnSim->setLevelTime(level, from_t, to_t);
  dustFluidSim->setLevelTime(level, from_t, to_t);
  //RK advance interior(including innner ghost cells) of level

#if USE_COSMOTRACE
  ray->preAdvance(hierarchy, ln);
#endif


  RKEvolveLevel(hierarchy, ln, from_t, to_t);


  level->getBoxLevel()->getMPI().Barrier();

  // recursively advancing children levels
  advanceLevel(hierarchy, ln+1, from_t, from_t + (to_t - from_t)/2.0);

  level->getBoxLevel()->getMPI().Barrier();
  
  advanceLevel(hierarchy, ln+1, from_t + (to_t - from_t)/2.0, to_t);

#if USE_COSMOTRACE
  ray->particleRedistribution(hierarchy, ln, true, step);
#endif
   
  // do some coarsening and
  // then update ghost cells through doing refinement if it has finer level
  if(ln < hierarchy->getNumberOfLevels() -1 )
  {
    level->getBoxLevel()->getMPI().Barrier();
    coarsen_schedules[ln]->coarsenData();

    level->getBoxLevel()->getMPI().Barrier();
    post_refine_schedules[ln]->fillData(to_t);
  }

#if USE_COSMOTRACE
  if(freeze_time_evolution == true && time_dependent_fields == true)
    bssnSim->set_time_dependent_fields(hierarchy, to_t);
#endif

  
  // copy _a to _p and set _p time to next timestamp

  math::HierarchyCellDataOpsReal<real_t> hcellmath(hierarchy,ln,ln);

  bssnSim->set_norm(level);
  
  bssnSim->copyAToP(hcellmath);
  dustFluidSim->copyAToP(hcellmath);
  
  bssnSim->setLevelTime(level, to_t, to_t);
  dustFluidSim->setLevelTime(level, to_t, to_t);
  
  level->getBoxLevel()->getMPI().Barrier();
}


void DustFluidSim::resetHierarchyConfiguration(
  /*! New hierarchy */
  const std::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
  /*! Coarsest level */ int coarsest_level,
  /*! Finest level */ int finest_level)
{
  pre_refine_schedules.resize(finest_level + 1);
  post_refine_schedules.resize(finest_level + 1);
  coarsen_schedules.resize(finest_level + 1);

  xfer::RefineAlgorithm pre_refiner, post_refiner;
  xfer::CoarsenAlgorithm coarsener(dim);  
  
  bssnSim->registerRKRefiner(pre_refiner, space_refine_op);
  bssnSim->registerCoarsenActive(coarsener,space_coarsen_op);
  bssnSim->registerRKRefinerActive(post_refiner, space_refine_op);

  dustFluidSim->registerRKRefinerActive(post_refiner, space_refine_op);
  dustFluidSim->registerRKRefiner(pre_refiner, space_refine_op);
  dustFluidSim->registerCoarsenActive(coarsener,space_coarsen_op);

  
  for(int ln = 0; ln <= finest_level; ln++)
  {
    const std::shared_ptr<hier::PatchLevel> level(
      new_hierarchy->getPatchLevel(ln));

    // reset pre refine refine schedule
    if(ln == 0)
    {
      pre_refine_schedules[ln] = pre_refiner.createSchedule(level, NULL);
    }
    else
    {
      pre_refine_schedules[ln] = pre_refiner.createSchedule(
        //border_fill_pattern,
        level,
        //level,
        ln - 1,
        new_hierarchy,
        (cosmoPS->is_time_dependent)?NULL:cosmoPS);
    }

    // reset coarse and post_refine schedule
    if(ln < finest_level)
    {
      coarsen_schedules[ln] = coarsener.createSchedule(level, new_hierarchy->getPatchLevel(ln+1));
      post_refine_schedules[ln] = post_refiner.createSchedule(level, NULL);
      
    }
  }
  
  return;
}
void DustFluidSim::putToRestart(
    const std::shared_ptr<tbox::Database>& restart_db) const
{
  restart_db->putDouble("cur_t", cur_t);
  restart_db->putInteger("step", step);
  return;
}

void DustFluidSim::getFromRestart()
{
  std::shared_ptr<tbox::Database> root_db(
    tbox::RestartManager::getManager()->getRootDatabase());

  if (!root_db->isDatabase(simulation_type)) {
    TBOX_ERROR("Restart database corresponding to "
               << simulation_type << " not found in restart file" << std::endl);
  }

  std::shared_ptr<tbox::Database> db(root_db->getDatabase(simulation_type));
  
  cur_t = db->getDouble("cur_t");

  starting_t = cur_t;

  step = db->getInteger("step");

  starting_step = step;

}

}
