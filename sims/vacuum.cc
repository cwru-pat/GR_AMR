#include "vacuum.h"
#include "../cosmo_includes.h"
#include "vacuum_macros.h"
#include "../utils/math.h"
#include "SAMRAI/xfer/PatchLevelBorderFillPattern.h"

using namespace SAMRAI;

namespace cosmo
{

/**
 * @brief Constructing VacuumSim object
 * 
 * @param hierarchy 
 * @param dimenstion
 * @param input database, which includes ALL sub databases
 * @param IO stream
 * @param name of simulation type
 * @param name of output file for VisIt
 */

VacuumSim::VacuumSim(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  const tbox::Dimension& dim_in,
  std::shared_ptr<tbox::InputDatabase>& input_db_in,
  std::ostream* l_stream_in,
  std::string simulation_type_in,
  std::string vis_filename_in):CosmoSim(
    hierarchy,
    dim_in, input_db_in, l_stream_in, simulation_type_in, vis_filename_in),
  cosmo_vacuum_db(input_db_in->getDatabase("VacuumSim"))
{

  t_init->start();

  std::string bd_type = cosmo_vacuum_db->getString("boundary_type");
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  
  gradient_indicator_idx =
    variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable(gradient_indicator), variable_db->getContext("ACTIVE"));

  
  // setting boundary type
  if(bd_type == "sommerfield")
  {
    cosmoPS = new SommerfieldBD(dim, bd_type);
  }
  else if(bd_type == "periodic")
  {
    cosmoPS = new periodicBD(dim, bd_type);
  }
  else
    TBOX_ERROR("Unsupported boundary type!\n");


  cosmo_vacuum_db = input_db->getDatabase("VacuumSim");
    

  tbox::pout<<"Running 'vacuum' type simulation.\n";

  // adding all fields to a list
  bssnSim->addFieldsToList(variable_id_list);

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

VacuumSim::~VacuumSim() {
}


void VacuumSim::getFromRestart()
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
void VacuumSim::init()
{
  
}

/**
 * @brief      Set vacuum initial conditions
 *
 */
void VacuumSim::setICs(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  tbox::plog<<"Setting initial conditions (ICs).";

  TBOX_ASSERT(gridding_algorithm);
  
  gridding_algorithm->printClassData(tbox::plog);

  bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();

  if(is_from_restart)
  {
    hierarchy->initializeHierarchy();
  }
 
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
  while(!is_from_restart && hierarchy->getNumberOfLevels() < hierarchy->getMaxNumberOfLevels())
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
  
  tbox::plog<<"Finished setting ICs. with hierarchy has "
            <<hierarchy->getNumberOfLevels()<<" levels\n";
}

void VacuumSim::initVacuumStep(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  bssnSim->stepInit(hierarchy);
}

/**
 * @brief      initilize newly created level
 *             set value directly if it's possible and return true
 *             do nothing when it's not possible and return false
 *
 */
bool VacuumSim::initLevel(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln)
{
  std::string ic_type = cosmo_vacuum_db->getString("ic_type");

  std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

  math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);

  // zero all fields
  bssnSim->clearSrc(hierarchy, ln);
  bssnSim->clearField(hierarchy, ln);
  bssnSim->clearGen1(hierarchy, ln);

  
  if(ic_type == "static_blackhole")
  {
    if(!USE_BSSN_SHIFT)
      TBOX_ERROR("Must enable shift for blackhole simulation!\n");
    bssn_ic_static_blackhole(hierarchy,ln);
    return true;
  }
  else if(ic_type == "kerr_blackhole")
  {
    if(!USE_BSSN_SHIFT)
      TBOX_ERROR("Must enable shift for blackhole simulation!\n");
    bssn_ic_kerr_blackhole(hierarchy,ln);
    return true;
  }
  else if(ic_type == "kerr_BHL_CTT")
  {
    if(ln > 0) return false;
    double M = cosmo_vacuum_db->getDoubleWithDefault("M", 1);
    double a = cosmo_vacuum_db->getDoubleWithDefault("spin", 0.5);
    double K_c = cosmo_vacuum_db->getDoubleWithDefault("K_c", 1);
    double relaxation_tolerance = cosmo_vacuum_db->getDoubleWithDefault("relaxation_tolerance", 1e-8);
    int num_vcycles = cosmo_vacuum_db->getIntegerWithDefault("num_vcycles", 4);
    int max_depth = cosmo_vacuum_db->getIntegerWithDefault("max_depth", 4);

    double l = cosmo_vacuum_db->getDoubleWithDefault("l", 1);
    double sigma = cosmo_vacuum_db->getDoubleWithDefault("sigma", 3.5);

    if(!USE_BSSN_SHIFT)
      TBOX_ERROR("Must enable shift for blackhole simulation!\n");
    
    bssn_ic_kerr_BHL_CTT(hierarchy,ln, M, a, K_c, relaxation_tolerance
                         , num_vcycles, max_depth, l, sigma);
    return true;
  }
  else if(ic_type == "static_BHL_CTT")
  {
    if(ln > 0) return false;
    double M = cosmo_vacuum_db->getDoubleWithDefault("M", 1);
    double K_c = cosmo_vacuum_db->getDoubleWithDefault("K_c", 1);
    double relaxation_tolerance = cosmo_vacuum_db->getDoubleWithDefault("relaxation_tolerance", 1e-8);
    int num_vcycles = cosmo_vacuum_db->getIntegerWithDefault("num_vcycles", 4);

    if(!USE_BSSN_SHIFT)
      TBOX_ERROR("Must enable shift for blackhole simulation!\n");
    
    bssn_ic_static_BHL_CTT(hierarchy,ln, M, 0, K_c, relaxation_tolerance, num_vcycles);
    return true;
  }
  else if(ic_type == "ds_blackhole")
  {
    if(!USE_BSSN_SHIFT)
      TBOX_ERROR("Must enable shift for blackhole simulation!\n");
    bssn_ic_ds_blackhole(hierarchy,ln,bssnSim);
    return true;
  }
  else if(ic_type == "awa_stability")
  {
    bssn_ic_awa_stability(hierarchy,ln,1e-10);
    return true;
  }
  else if(ic_type == "awa_linear_wave")
  {
    
    bssn_ic_awa_linear_wave(hierarchy,ln,1e-8,1);
    return true;
  }
  else if(ic_type == "awa_gauge_wave")
  {
    bssn_ic_awa_gauge_wave(hierarchy,ln,1);
    return true;
  }
  else if(ic_type == "awa_shifted_gauge_wave")
  {
    bssn_ic_awa_shifted_gauge_wave(hierarchy,ln,1);
    return true;
  }

  else
    TBOX_ERROR("Undefined IC type!\n");

  return false;
}

/**
 * @brief      compute weight for each grid, corresponds grid volumn, 
 *             equals to 0 if it's covered by finer level
 *
 */  
void VacuumSim::computeVectorWeights(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  TBOX_ASSERT(hierarchy);
  TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, *hierarchy);

  int weight_id = weight_idx;
  int coarsest_ln = 0;
  int finest_ln = hierarchy->getFinestLevelNumber();


  int ln;
  for (ln = finest_ln; ln >= coarsest_ln; --ln) {

    /*
     * On every level, first assign cell volume to vector weight.
     */

    std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    for (hier::PatchLevel::iterator p(level->begin());
         p != level->end(); ++p) {
      const std::shared_ptr<hier::Patch>& patch = *p;
      std::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      TBOX_ASSERT(patch_geometry);

      const double* dx = patch_geometry->getDx();
      double cell_vol = dx[0];
      if (dim > tbox::Dimension(1)) {
        cell_vol *= dx[1];
      }

      if (dim > tbox::Dimension(2)) {
        cell_vol *= dx[2];
      }

      std::shared_ptr<pdat::CellData<double> > w(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_id)));
      TBOX_ASSERT(w);
      w->fillAll(cell_vol);
    }

    /*
     * On all but the finest level, assign 0 to vector
     * weight to cells covered by finer cells.
     */

    if (ln < finest_ln) {

      /*
       * First get the boxes that describe index space of the next finer
       * level and coarsen them to describe corresponding index space
       * at this level.
       */

      std::shared_ptr<hier::PatchLevel> next_finer_level(
        hierarchy->getPatchLevel(ln + 1));
      hier::BoxContainer coarsened_boxes = next_finer_level->getBoxes();
      hier::IntVector coarsen_ratio(next_finer_level->getRatioToLevelZero());
      coarsen_ratio /= level->getRatioToLevelZero();
      coarsened_boxes.coarsen(coarsen_ratio);

      /*
       * Then set vector weight to 0 wherever there is
       * a nonempty intersection with the next finer level.
       * Note that all assignments are local.
       */

      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {

        const std::shared_ptr<hier::Patch>& patch = *p;
        for (hier::BoxContainer::iterator i = coarsened_boxes.begin();
             i != coarsened_boxes.end(); ++i) {

          hier::Box intersection = *i * (patch->getBox());
          if (!intersection.empty()) {
            std::shared_ptr<pdat::CellData<double> > w(
              SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(weight_id)));
            TBOX_ASSERT(w);
            w->fillAll(0.0, intersection);

          }  // assignment only in non-empty intersection
        }  // loop over coarsened boxes from finer level
      }  // loop over patches in level
    }  // all levels except finest
  }  // loop over levels
  
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
void VacuumSim::initializeLevelData(
   /*! Hierarchy to initialize */
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   /*! Level to initialize */
   const int ln,
   const double init_data_time,
   const bool can_be_refined,
   /*! Whether level is being introduced for the first time */
   const bool initial_time,
   /*! Level to copy data from */
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
     level->allocatePatchData(weight_idx);
     // if(use_AHFinder)
     //   horizon->alloc(patch_hierarchy, ln);
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
   // if(use_AHFinder)
   //   horizon->clear(patch_hierarchy, ln);
   /*
    * Refine solution data from coarser level and, if provided, old level.
    */
   
   xfer::RefineAlgorithm refiner;

   std::shared_ptr<hier::RefineOperator> accurate_refine_op =
     space_refine_op;
     
   TBOX_ASSERT(accurate_refine_op);

   //registering refine variables
   bssnSim->registerRKRefinerActive(refiner, accurate_refine_op);
                             
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

   // if(use_AHFinder)
   //   horizon->copyAToP(hcellmath);
   
   level->getBoxLevel()->getMPI().Barrier();
   /* Set vector weight. */
   computeVectorWeights(hierarchy);
}

/**
 * @brief tag the grids that need to be refined
 *
 */  
void VacuumSim::applyGradientDetector(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy_,
   const int ln,
   const double error_data_time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation)
{
  // do not tag new grid when restarting
  if(tbox::RestartManager::getManager()->isFromRestart() && step == starting_step)
    return;
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


void VacuumSim::outputVacuumStep(
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
  
  bssnSim->output_L2_H_constaint(
    hierarchy, weight_idx, cosmoPS, 0.5);
  //  bssnSim->output_max_H_constaint(hierarchy, weight_idx);
 
  cosmo_io->registerVariablesWithPlotter(*visit_writer, step);
  cosmo_io->dumpData(hierarchy, *visit_writer, step, cur_t);

  cosmo_statistic->output_expansion_info(
      hierarchy,
      bssnSim, weight_idx, step, cur_t, max_horizon_radius);

  cosmo_statistic->output_conformal_avg(
    hierarchy,
    bssnSim, weight_idx, step, cur_t, max_horizon_radius);

}

/**
 * @brief advance hierarchy from time "from_t" to "to_t"
 * 
 * @param hierarchy
 * @param starting time
 * @param ending time
 */  
void VacuumSim::runVacuumStep(
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
double VacuumSim::getDt(
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
void VacuumSim::runStep(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  //  std::cout<<cur_t<<" ";
  runCommonStepTasks(hierarchy);

  double dt = getDt(hierarchy);
  initVacuumStep(hierarchy);
  
  outputVacuumStep(hierarchy);
  runVacuumStep(hierarchy, cur_t, cur_t + dt);
  cur_t += dt;
}


void VacuumSim::addBSSNExtras(
  const std::shared_ptr<hier::PatchLevel> & level)
{
  return;
}

void VacuumSim::addBSSNExtras(
  const std::shared_ptr<hier::Patch> & patch)
{
  return;
}

/**
 * @brief RK evolve level
 * 
 * @param hierarchy
 * @param level index
 * @param starting time
 * @param ending time
 */  
void VacuumSim::RKEvolveLevel(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln,
  double from_t,
  double to_t)
{
  const std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  const std::shared_ptr<hier::PatchLevel> coarser_level(
    ((ln>0)?(hierarchy->getPatchLevel(ln-1)):NULL));
  
  

  bssnSim->prepareForK1(coarser_level, to_t);
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    //Evolve inner grids
    bssnSim->RKEvolvePatch(patch, to_t - from_t);
    // Evolve physical boundary
    // would not do anything if boundary is time independent

    bssnSim->RKEvolvePatchBD(patch, to_t - from_t);
  }

  
  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
  pre_refine_schedules[ln]->fillData(to_t);

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
    bssnSim->K1FinalizePatch(patch);    
    addBSSNExtras(patch);
  }
  bssnSim->set_norm(level);
  
  
  /**************Starting K2 *********************************/
  bssnSim->prepareForK2(coarser_level, to_t);
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    //Evolve inner grids
    bssnSim->RKEvolvePatch(patch, to_t - from_t);
    //Evolve physical boundary
    // would not do anything if boundary is time dependent
    bssnSim->RKEvolvePatchBD(patch, to_t - from_t);  
  }

  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
  pre_refine_schedules[ln]->fillData(to_t);

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
    bssnSim->K2FinalizePatch(patch);
    addBSSNExtras(patch);
  }

  bssnSim->set_norm(level);

  /**************Starting K3 *********************************/

  bssnSim->prepareForK3(coarser_level, to_t);
    

  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    bssnSim->RKEvolvePatch(patch, to_t - from_t);
    //Evolve physical boundary
    // would not do anything if boundary is time dependent
    bssnSim->RKEvolvePatchBD(patch, to_t - from_t);  
  }

  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
  pre_refine_schedules[ln]->fillData(to_t);

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
    bssnSim->K3FinalizePatch(patch);
    addBSSNExtras(patch);
  }

  bssnSim->set_norm(level);

  /**************Starting K4 *********************************/

  bssnSim->prepareForK4(coarser_level, to_t);

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    bssnSim->RKEvolvePatch(patch, to_t - from_t);
    //Evolve physical boundary
    // would not do anything if boundary is time dependent
    bssnSim->RKEvolvePatchBD(patch, to_t - from_t);  
  }

  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
  pre_refine_schedules[ln]->fillData(to_t);

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
    bssnSim->K4FinalizePatch(patch);
    addBSSNExtras(patch);
  }

  bssnSim->set_norm(level);

}

/**
 * @brief advance a single level by:
 *        1. RK evolve this level(including evolve the interior and interpolate the boundary)
 *        2. evolve its son levels recursively
 *        3. doing coarsen operation
 */ 
void VacuumSim::advanceLevel(
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
  //RK advance interior(including innner ghost cells) of level


  RKEvolveLevel(hierarchy, ln, from_t, to_t);


  level->getBoxLevel()->getMPI().Barrier();

  // recursively advancing children levels
  advanceLevel(hierarchy, ln+1, from_t, from_t + (to_t - from_t)/2.0);

  level->getBoxLevel()->getMPI().Barrier();
  
  advanceLevel(hierarchy, ln+1, from_t + (to_t - from_t)/2.0, to_t);
   
  // do some coarsening and
  // then update ghost cells through doing refinement if it has finer level
  if(ln < hierarchy->getNumberOfLevels() -1 )
  {
    level->getBoxLevel()->getMPI().Barrier();
    coarsen_schedules[ln]->coarsenData();

    level->getBoxLevel()->getMPI().Barrier();
    post_refine_schedules[ln]->fillData(to_t);
  }

  // copy _a to _p and set _p time to next timestamp

  math::HierarchyCellDataOpsReal<real_t> hcellmath(hierarchy,ln,ln);

  bssnSim->set_norm(level);
  
  bssnSim->copyAToP(hcellmath);

  bssnSim->setLevelTime(level, to_t, to_t);

  level->getBoxLevel()->getMPI().Barrier();
}


void VacuumSim::resetHierarchyConfiguration(
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
  
void VacuumSim::putToRestart(
    const std::shared_ptr<tbox::Database>& restart_db) const
{
  restart_db->putDouble("cur_t", cur_t);
  restart_db->putInteger("step", step);
  return;
}


} /* namespace cosmo */
