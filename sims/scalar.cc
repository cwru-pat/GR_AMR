#include "scalar.h"
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

ScalarSim::ScalarSim(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  const tbox::Dimension& dim_in,
  boost::shared_ptr<tbox::InputDatabase>& input_db_in,
  std::ostream* l_stream_in,
  std::string simulation_type_in,
  std::string vis_filename_in):CosmoSim(
    hierarchy,
    dim_in, input_db_in, l_stream_in, simulation_type_in, vis_filename_in),
  cosmo_scalar_db(input_db_in->getDatabase("ScalarSim"))
{
  stop_after_setting_init =
    cosmo_scalar_db->getBoolWithDefault("stop_after_setting_init", false);

  t_init->start();

  std::string bd_type = cosmo_scalar_db->getString("boundary_type");

  if(bd_type == "periodic")
  {
    cosmoPS = new periodicBD(dim, bd_type);
  }
  else
    TBOX_ERROR("Unsupported boundary type!\n");


  cosmo_scalar_db = input_db->getDatabase("ScalarSim");
    

  tbox::pout<<"Running 'scalar' type simulation.\n";

  // adding all fields to a list
  bssnSim->addFieldsToList(variable_id_list);

  scalarSim = new Scalar(
    hierarchy, dim, input_db->getDatabase("Scalar"), lstream,KO_damping_coefficient);

  //variable_id_list.push_back(staticSim->DIFFD_a_idx);

  scalarSim->addFieldsToList(variable_id_list);
  
  variable_id_list.push_back(weight_idx);

  hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);
  
  t_init->stop();  
}

ScalarSim::~ScalarSim() {
}

  
void ScalarSim::init()
{
  
}

/**
 * @brief      Set vacuum initial conditions
 *
 */
void ScalarSim::setICs(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  tbox::plog<<"Setting initial conditions (ICs).";

  TBOX_ASSERT(gridding_algorithm);
  
  gridding_algorithm->printClassData(tbox::plog);

  // set initial condition by calling function initializeLevelData()
  gridding_algorithm->makeCoarsestLevel(0.0);

  // regrid initial hierarchy if needed
  while(hierarchy->getNumberOfLevels() < hierarchy->getMaxNumberOfLevels())
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
      0.0);
    int post_level_num = hierarchy->getNumberOfLevels();
    // no new level is created
    if(post_level_num == pre_level_num) break;
  }
  
  tbox::plog<<"Finished setting ICs. with hierarchy has "
            <<hierarchy->getNumberOfLevels()<<" levels\n";
}

void ScalarSim::initScalarStep(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  bssnSim->stepInit(hierarchy);
  bssnSim->clearSrc(hierarchy);
  scalarSim->addBSSNSrc(bssnSim, hierarchy);
  if(stop_after_setting_init)
    TBOX_ERROR("Stop after initializing hierarchy as demanded\n");
}

/**
 * @brief      initilize newly created level
 *             set value directly if it's possible and return true
 *             do nothing when it's not possible and return false
 *
 */
bool ScalarSim::initLevel(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln)
{
  std::string ic_type = cosmo_scalar_db->getString("ic_type");

  boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

  math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);

  // zero all fields
  bssnSim->clearField(hierarchy, ln);
  bssnSim->clearSrc(hierarchy, ln);
  bssnSim->clearGen1(hierarchy, ln);

  scalarSim->clear(hierarchy, ln);
  
  
  if(ic_type == "semianalytic_test")
  {
    if(ln > 0) return false;
    scalar_ic_set_semianalytic_test(hierarchy, ln, bssnSim, scalarSim, input_db->getDatabase("Scalar"));
    return true;
  }
  else if(ic_type == "scalar_collapse")
  {
    // which means not initial data file exist
    scalar_ic_set_scalar_collapse(hierarchy, ln, bssnSim, scalarSim, input_db->getDatabase("Scalar"));
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
void ScalarSim::computeVectorWeights(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
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

    boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    for (hier::PatchLevel::iterator p(level->begin());
         p != level->end(); ++p) {
      const boost::shared_ptr<hier::Patch>& patch = *p;
      
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
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

      boost::shared_ptr<pdat::CellData<double> > w(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
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

      boost::shared_ptr<hier::PatchLevel> next_finer_level(
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

        const boost::shared_ptr<hier::Patch>& patch = *p;
        for (hier::BoxContainer::iterator i = coarsened_boxes.begin();
             i != coarsened_boxes.end(); ++i) {

          hier::Box intersection = *i * (patch->getBox());
          if (!intersection.empty()) {
            boost::shared_ptr<pdat::CellData<double> > w(
              BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
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

void ScalarSim::initializeLevelData(
   /*! Hierarchy to initialize */
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   /*! Level to initialize */
   const int ln,
   const double init_data_time,
   const bool can_be_refined,
   /*! Whether level is being introduced for the first time */
   const bool initial_time,
   /*! Level to copy data from */
   const boost::shared_ptr<hier::PatchLevel>& old_level,
   const bool allocate_data)
{
   NULL_USE(can_be_refined);
   NULL_USE(initial_time);

   boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy(hierarchy);

   /*
    * Reference the level object with the given index from the hierarchy.
    */
   boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
   

   math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);
      
   
   if (allocate_data)
   {
     bssnSim->allocField(patch_hierarchy, ln);
     bssnSim->allocSrc(patch_hierarchy, ln);
     bssnSim->allocGen1(patch_hierarchy, ln);
     scalarSim->alloc(patch_hierarchy, ln);
     level->allocatePatchData(weight_idx);
   }

   // marks whether we have solved initial value for certain level,
   // if so, there is no need to interpolate from coarser level
   bool has_initial = false;

   //at beginning, initialize new level
   if(init_data_time < EPS)
   {
     has_initial = initLevel(patch_hierarchy, ln);
   }
   bssnSim->clearSrc(patch_hierarchy, ln);
   bssnSim->clearGen1(patch_hierarchy, ln);
   /*
    * Refine solution data from coarser level and, if provided, old level.
    */
   
   xfer::RefineAlgorithm refiner;

   boost::shared_ptr<hier::RefineOperator> accurate_refine_op =
     space_refine_op;
     
   TBOX_ASSERT(accurate_refine_op);

   //registering refine variables
   //BSSN_APPLY_TO_FIELDS_ARGS(DUS_REGISTER_SPACE_REFINE_A,refiner,accurate_refine_op);

   bssnSim->registerRKRefinerActive(refiner, accurate_refine_op);
   scalarSim->registerRKRefinerActive(refiner, accurate_refine_op);
   
   boost::shared_ptr<xfer::RefineSchedule> refine_schedule;

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
   scalarSim->copyAToP(hcellmath);
   
   level->getBoxLevel()->getMPI().Barrier();
   /* Set vector weight. */
   computeVectorWeights(hierarchy);
}

/**
 * @brief tag the grids that need to be refined
 *
 */  
void ScalarSim::applyGradientDetector(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy_,
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
   boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
       hierarchy.getGridGeometry()));
   double max_der_norm = 0;
   hier::PatchLevel& level =
      (hier::PatchLevel &) * hierarchy.getPatchLevel(ln);
   size_t ntag = 0, ntotal = 0;
   //double maxestimate = 0;
   for (hier::PatchLevel::iterator pi(level.begin());
        pi != level.end(); ++pi)
   {
      hier::Patch& patch = **pi;

      const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch.getPatchGeometry()));

      boost::shared_ptr<hier::PatchData> tag_data(
         patch.getPatchData(tag_index));
      ntotal += patch.getBox().numberCells().getProduct();
      if (!tag_data)
      {
         TBOX_ERROR(
            "Data index " << tag_index << " does not exist for patch.\n");
      }
      boost::shared_ptr<pdat::CellData<int> > tag_cell_data_(
         BOOST_CAST<pdat::CellData<int>, hier::PatchData>(tag_data));
      TBOX_ASSERT(tag_cell_data_);
      
      boost::shared_ptr<pdat::CellData<double>> K_data(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
          patch.getPatchData(bssnSim->DIFFchi_a_idx)));


      arr_t K = pdat::ArrayDataAccess::access<DIM, real_t>(
        K_data->getArrayData());
      
      if (!K_data) {
         TBOX_ERROR("Data index " << bssnSim->DIFFchi_p_idx
                                  << " does not exist for patch.\n");
      }
      pdat::CellData<idx_t>& tag_cell_data = *tag_cell_data_;
          
      tag_cell_data.fill(0);
      
      hier::Box::iterator iend(patch.getBox().end());

      for (hier::Box::iterator i(patch.getBox().begin()); i != iend; ++i)
      {
         const pdat::CellIndex cell_index(*i);
         max_der_norm = tbox::MathUtilities<double>::Max(
           max_der_norm,
           derivative_norm(
             cell_index(0),
             cell_index(1),
             cell_index(2),
             K));
         if(derivative_norm(
              cell_index(0),
              cell_index(1),
              cell_index(2),
              K) > adaption_threshold)
         {
          
           tag_cell_data(cell_index) = 1;
           ++ntag;
         }
       
      }

   }
   const tbox::SAMRAI_MPI& mpi(hierarchy.getMPI());
   if (mpi.getSize() > 1)
   {
     mpi.AllReduce(&max_der_norm, 1, MPI_MAX);
   }

   tbox::plog << "Adaption threshold is " << adaption_threshold << "\n";
   tbox::plog << "Number of cells tagged on level " << ln << " is "
              << ntag << "/" << ntotal << "\n";
   tbox::plog << "Max norm is " << max_der_norm << "\n";
}
  
void ScalarSim::outputScalarStep(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  boost::shared_ptr<appu::VisItDataWriter> visit_writer(
    new appu::VisItDataWriter(
    dim, "VisIt Writer", vis_filename + ".visit"));

  tbox::pout<<"step: "<<step<<"/"<<num_steps<<"\n";
  
  bssnSim->output_L2_H_constaint(hierarchy, weight_idx);
 
  cosmo_io->registerVariablesWithPlotter(*visit_writer, step);
  cosmo_io->dumpData(hierarchy, *visit_writer, step, cur_t);
  
}

/**
 * @brief advance hierarchy from time "from_t" to "to_t"
 * 
 * @param hierarchy
 * @param starting time
 * @param ending time
 */  
void ScalarSim::runScalarStep(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
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
double ScalarSim::getDt(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
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
void ScalarSim::runStep(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  runCommonStepTasks(hierarchy);
  double dt = getDt(hierarchy);
  initScalarStep(hierarchy);
  
  outputScalarStep(hierarchy);

  runScalarStep(hierarchy, cur_t, cur_t + dt);
  cur_t += dt;
}


void ScalarSim::addBSSNExtras(
  const boost::shared_ptr<hier::PatchLevel> & level)
{
  return;
}

void ScalarSim::addBSSNExtras(
  const boost::shared_ptr<hier::Patch> & patch)
{
  return;
}

void ScalarSim::RKEvolve(
  const boost::shared_ptr<hier::Patch> & patch, real_t dt)
{
  bssnSim->initPData(patch);
  bssnSim->initMDA(patch);
  scalarSim->initPData(patch);
  scalarSim->initMDA(patch);
  const hier::Box& box = patch->getBox();
  
  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  BSSNData bd = {0};
  
  const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(  
    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  //initialize dx for each patch
  const real_t * dx = &(patch_geom->getDx())[0];


  
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        bssnSim->RKEvolvePt(i, j, k, bd, dx, dt);
        scalarSim->RKEvolvePt(i, j, k, bd, dx, dt);
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
void ScalarSim::RKEvolveLevel(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln,
  double from_t,
  double to_t)
{
  const boost::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  const boost::shared_ptr<hier::PatchLevel> coarser_level(
    ((ln>0)?(hierarchy->getPatchLevel(ln-1)):NULL));
  
  
  xfer::RefineAlgorithm refiner;
  boost::shared_ptr<xfer::RefineSchedule> refine_schedule;


  
  bssnSim->registerRKRefiner(refiner, space_refine_op);
  scalarSim->registerRKRefiner(refiner, space_refine_op);
  
  bssnSim->prepareForK1(coarser_level, to_t);
  scalarSim->prepareForK1(coarser_level, to_t);

  
  //if not the coarsest level, should 
  if(coarser_level!=NULL)
  {
    boost::shared_ptr<xfer::PatchLevelBorderFillPattern> border_fill_pattern (
      new xfer::PatchLevelBorderFillPattern());
    
    refine_schedule = refiner.createSchedule(
      //border_fill_pattern,
      level,
      //level,
      coarser_level->getLevelNumber(),
      hierarchy,
      (cosmoPS->is_time_dependent)?NULL:cosmoPS);
  }
  else
    refine_schedule = refiner.createSchedule(level, NULL);
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;

    //Evolve inner grids
    RKEvolve(patch, to_t - from_t);
  }

  
  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
  refine_schedule->fillData(to_t);

  bssnSim->clearSrc(hierarchy, ln);
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    bssnSim->K1FinalizePatch(patch);
    scalarSim->K1FinalizePatch(patch);
    scalarSim->addBSSNSrc(bssnSim,patch);
  }
  bssnSim->set_norm(level);
  
  /**************Starting K2 *********************************/
  bssnSim->prepareForK2(coarser_level, to_t);
  scalarSim->prepareForK2(coarser_level, to_t);
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;

    //Evolve inner grids
    RKEvolve(patch, to_t - from_t);
  }

  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
  refine_schedule->fillData(to_t);

  bssnSim->clearSrc(hierarchy, ln);
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    bssnSim->K2FinalizePatch(patch);
    scalarSim->K2FinalizePatch(patch);
    scalarSim->addBSSNSrc(bssnSim, patch);
  }

  bssnSim->set_norm(level);
  
  /**************Starting K3 *********************************/

  bssnSim->prepareForK3(coarser_level, to_t);
  scalarSim->prepareForK3(coarser_level, to_t);    

  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    RKEvolve(patch, to_t - from_t);

  }

  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
  refine_schedule->fillData(to_t);

  bssnSim->clearSrc(hierarchy, ln);
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    bssnSim->K3FinalizePatch(patch);
    scalarSim->K3FinalizePatch(patch);
    scalarSim->addBSSNSrc(bssnSim,patch);
  }

  bssnSim->set_norm(level);
  
  /**************Starting K4 *********************************/

  bssnSim->prepareForK4(coarser_level, to_t);
  scalarSim->prepareForK4(coarser_level, to_t);

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    RKEvolve(patch, to_t - from_t);

  }

  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
  refine_schedule->fillData(to_t);

  bssnSim->clearSrc(hierarchy, ln);
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    bssnSim->K4FinalizePatch(patch);
    scalarSim->K4FinalizePatch(patch);
  }
  bssnSim->set_norm(level);
}

/**
 * @brief advance a single level by:
 *        1. RK evolve this level(including evolve the interior and interpolate the boundary)
 *        2. evolve its son levels recursively
 *        3. doing coarsen operation
 */ 
void ScalarSim::advanceLevel(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  int ln,
  double from_t,
  double to_t)
{
  if( ln >= hierarchy->getNumberOfLevels())
    return;
  
  //double dt = to_t - from_t;
  const boost::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  //updating extra fields before advancing any level
  addBSSNExtras(level);

  bssnSim->setLevelTime(level, from_t, to_t);
  scalarSim->setLevelTime(level, from_t, to_t);
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
    xfer::CoarsenAlgorithm coarsener(dim);
  

    boost::shared_ptr<xfer::CoarsenSchedule> coarsen_schedule;

    bssnSim->registerCoarsenActive(coarsener,space_coarsen_op);
    scalarSim->registerCoarsenActive(coarsener,space_coarsen_op);
      
    coarsen_schedule = coarsener.createSchedule(level, hierarchy->getPatchLevel(ln+1));
    level->getBoxLevel()->getMPI().Barrier();
    coarsen_schedule->coarsenData();

    xfer::RefineAlgorithm post_refiner;

    boost::shared_ptr<xfer::RefineSchedule> refine_schedule;

    

    bssnSim->registerRKRefinerActive(post_refiner, space_refine_op);
    scalarSim->registerRKRefinerActive(post_refiner, space_refine_op);
    
    refine_schedule = post_refiner.createSchedule(level, NULL);

    level->getBoxLevel()->getMPI().Barrier();
    refine_schedule->fillData(to_t);
  }

  // copy _a to _p and set _p time to next timestamp

  math::HierarchyCellDataOpsReal<real_t> hcellmath(hierarchy,ln,ln);

  bssnSim->set_norm(level);
  
  bssnSim->copyAToP(hcellmath);
  scalarSim->copyAToP(hcellmath);

  bssnSim->setLevelTime(level, to_t, to_t);
  scalarSim->setLevelTime(level, to_t, to_t);

  level->getBoxLevel()->getMPI().Barrier();
}


void ScalarSim::resetHierarchyConfiguration(
  /*! New hierarchy */
  const boost::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
  /*! Coarsest level */ int coarsest_level,
  /*! Finest level */ int finest_level)
{
  return;
}
  

}
