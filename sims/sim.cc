#include "sim.h"
#include "../cosmo_includes.h"
#include "../components/bssn/bssn.h"
#include "../utils/math.h"
#include "SAMRAI/xfer/PatchLevelBorderFillPattern.h"
#include "SAMRAI/tbox/MathUtilities.h"

using namespace SAMRAI;


namespace cosmo
{
/*
 * Claming timers 
 */  
boost::shared_ptr<tbox::Timer> CosmoSim::t_loop;
boost::shared_ptr<tbox::Timer> CosmoSim::t_init;
boost::shared_ptr<tbox::Timer> CosmoSim::t_RK_steps;

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
CosmoSim::CosmoSim(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  const tbox::Dimension& dim_in,
  boost::shared_ptr<tbox::InputDatabase>& input_db_in,
  std::ostream* l_stream_in = 0,
  std::string simulation_type_in = std::string(),
  std::string vis_filename_in = std::string()):
  input_db(input_db_in),
  cosmo_sim_db(input_db->getDatabase("CosmoSim")),
  lstream(l_stream_in),
  dim(dim_in),
  step(0),
  starting_step(0),
  num_steps(cosmo_sim_db->getInteger("steps")),
  simulation_type(simulation_type_in),
  comments(cosmo_sim_db->getStringWithDefault("comments","")),
  do_plot(cosmo_sim_db->getBoolWithDefault("do_plot", false)),
  dt_frac(cosmo_sim_db->getDoubleWithDefault("dt_frac", 0.1)),
  vis_filename(vis_filename_in),
  cur_t(0),
  starting_t(0),
  weight(new pdat::CellVariable<real_t>(dim, "weight", 1)),
  refine_scratch(new pdat::CellVariable<real_t>(dim, "refine_scratch", 1)),
  weight_idx(0),
  regridding_interval(cosmo_sim_db->getInteger("regridding_interval")),
  KO_damping_coefficient(cosmo_sim_db->getDoubleWithDefault("KO_damping_coefficient",0)),
  adaption_threshold(cosmo_sim_db->getDoubleWithDefault("adaption_threshold", 1)),
  refine_op_type(cosmo_sim_db->getStringWithDefault("refine_op_type", "LINEAR_REFINE")),
  coarsen_op_type(cosmo_sim_db->getStringWithDefault("coarsen_op_type", "CONSERVATIVE_COARSEN")),
  use_AHFinder(cosmo_sim_db->getBoolWithDefault("use_AHFinder", false)),
  AHFinder_iter_limit(cosmo_sim_db->getIntegerWithDefault("AHFinder_iter_limit", 100)),
  AHFinder_dt_frac(cosmo_sim_db->getDoubleWithDefault("AHFinder_dt_frac", 0.1)),
  surface_move_shreshold(cosmo_sim_db->getDoubleWithDefault("surface_move_shreshold", 1e-9)),
  save_interval(cosmo_sim_db->getIntegerWithDefault("save_interval", std::numeric_limits<int>::max())),
  use_anguler_momentum_finder(cosmo_sim_db->getBoolWithDefault("use_anguler_momentum_finder", false))
{
  t_loop = tbox::TimerManager::getManager()->
    getTimer("loop");
  t_init = tbox::TimerManager::getManager()->
    getTimer("init");
  t_RK_steps = tbox::TimerManager::getManager()->
    getTimer("RK_steps");

  // initialzing BSSN object
  bssnSim = new BSSN(
    hierarchy, dim,input_db->getDatabase("BSSN"), lstream,KO_damping_coefficient);

  // initializing IO object
  cosmo_io = new CosmoIO(dim, input_db->getDatabase("IO"), lstream);

  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

  // get context ACTIVE to initilize these two extra fields 
  boost::shared_ptr<hier::VariableContext> context_active(
    variable_db->getContext("ACTIVE"));

  // weight component stores volume for each cell
  weight_idx = variable_db->registerVariableAndContext( 
      weight, 
      context_active,
      hier::IntVector(dim, STENCIL_ORDER));

  horizon = new Horizon(
    hierarchy, dim,input_db->getDatabase("Horizon"), lstream, weight_idx);

  
  // scractch component for refinement, currently not in use
  refine_scratch_idx = variable_db->registerVariableAndContext( 
      refine_scratch, 
      context_active,
      hier::IntVector(dim, STENCIL_ORDER));

  AHFinder_interval = cosmo_sim_db->getIntegerWithDefault("AHFinder_interval" , 0);

  
}

CosmoSim::~CosmoSim()
{

}
  
/**
 * @brief  set regriding algorithm
 */
void CosmoSim::setGriddingAlgs(
  boost::shared_ptr<mesh::GriddingAlgorithm>& gridding_algorithm_in)
{
  gridding_algorithm = gridding_algorithm_in;
}
/**
 * @brief  set refine and coarsen operators
 */
void CosmoSim::setRefineCoarsenOps(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{

  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));

  TBOX_ASSERT(grid_geometry_);

  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  space_refine_op =
    grid_geometry.
    lookupRefineOperator(bssnSim->DIFFchi, refine_op_type);
  
  TBOX_ASSERT(space_refine_op);

  space_coarsen_op =
    grid_geometry.
    lookupCoarsenOperator(bssnSim->DIFFchi, coarsen_op_type);

  TBOX_ASSERT(space_coarsen_op);
}

/**
 * @brief      Run the simulation.
 */
void CosmoSim::run(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  tbox::plog<<"Running simulation...";

  t_loop->start();
  while(step <= num_steps)
  {
    runStep(hierarchy);
    step++;
    
  }
  t_loop->stop();

  tbox::plog<<"\nEnding simulation.";
}



/**
 * @brief regrid when necessary and detect NaNs.
 */
void CosmoSim::runCommonStepTasks(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
    // detecting NaNs
  isValid(hierarchy);
  // since all neccecery levels were built when setting IC
  // no need to regrid again at zero step
  if(step > starting_step && (step % regridding_interval == 0))
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
    tbox::plog << "Newly adapted hierarchy\n";
    hierarchy->recursivePrint(tbox::plog, "    ", 1);
  }

  if(step > starting_step && (step %save_interval == 0))
  {
    std::string restart_file_name = simulation_type + comments +".restart";
    tbox::RestartManager::getManager()->writeRestartFile(restart_file_name, step);
  }

  bool found_horizon = false;
  
  if( step > 0 &&
     use_AHFinder && (step % AHFinder_interval == 0))
  {

    // starting finding apparent horizon
    horizon->initSurface(hierarchy);

    for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
    {
      boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
   

      math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);

      xfer::RefineAlgorithm refiner;

      boost::shared_ptr<hier::RefineOperator> accurate_refine_op =
        space_refine_op;
     

      //registering refine variables
      horizon->registerRKRefinerActive(refiner, accurate_refine_op);
                             
      boost::shared_ptr<xfer::RefineSchedule> refine_schedule;

      level->getBoxLevel()->getMPI().Barrier();
      refine_schedule =
        refiner.createSchedule(level,
                               level,
                               NULL);
     
   
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
 
      horizon->copyAToP(hcellmath);
    }
    
    found_horizon = findHorizon(hierarchy);
  }
  if(found_horizon && use_anguler_momentum_finder)
  {
    horizon->findKilling(hierarchy, bssnSim);
  }
}

void CosmoSim::initHorizonStep(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  horizon->addNormVector(bssnSim, hierarchy);
}

  
void CosmoSim::RKEvolveHorizonLevel(
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

  
  horizon->registerRKRefiner(refiner, space_refine_op);
  
  horizon->prepareForK1(coarser_level, to_t);
  
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
      NULL);
  }
  else
    refine_schedule = refiner.createSchedule(level, NULL);

  
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    //Evolve inner grids
    horizon->RKEvolveHorizon(patch, bssnSim, to_t - from_t);

  }
  
  
  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
  refine_schedule->fillData(to_t);


  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    horizon->K1FinalizePatch(patch);
    horizon->updateBD(patch, to_t - from_t);

    horizon->addNormVector(bssnSim,patch);
  }


  
  /**************Starting K2 *********************************/
  horizon->prepareForK2(coarser_level, to_t);

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    
    //Evolve inner grids
    horizon->RKEvolveHorizon(patch, bssnSim, to_t - from_t);
    
    
  }

  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
  refine_schedule->fillData(to_t);

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    horizon->K2FinalizePatch(patch);
    horizon->updateBD(patch, to_t - from_t);
        
    horizon->addNormVector(bssnSim,patch);
  }

  
  
  /**************Starting K3 *********************************/

  horizon->prepareForK3(coarser_level, to_t);

  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    horizon->RKEvolveHorizon(patch, bssnSim, to_t - from_t);

  }

  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
  refine_schedule->fillData(to_t);

  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    horizon->K3FinalizePatch(patch);
    horizon->updateBD(patch, to_t - from_t);
        
    horizon->addNormVector(bssnSim,patch);
  }


  
  /**************Starting K4 *********************************/

  horizon->prepareForK4(coarser_level, to_t);

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    horizon->RKEvolveHorizon(patch, bssnSim, to_t - from_t);
  }

  // fill ghost cells 
  level->getBoxLevel()->getMPI().Barrier();
  refine_schedule->fillData(to_t);

  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    horizon->K4FinalizePatch(patch);
    horizon->updateBD(patch, to_t - from_t);
        
    horizon->addNormVector(bssnSim,patch);
  }
  
}


bool CosmoSim::advanceHorizonLevel(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  int ln,
  double from_t,
  double to_t)
{
  if( ln >= hierarchy->getNumberOfLevels())
    return false;

    
  bool surface_is_moved = false;
  
  //double dt = to_t - from_t;
  const boost::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  horizon->setLevelTime(level, from_t, to_t);


  RKEvolveHorizonLevel(hierarchy, ln, from_t, to_t);

  level->getBoxLevel()->getMPI().Barrier();

  // recursively advancing children levels
  if(advanceHorizonLevel(hierarchy, ln+1, from_t, from_t + (to_t - from_t)/2.0))
    surface_is_moved = true;

  level->getBoxLevel()->getMPI().Barrier();
  
  if(advanceHorizonLevel(hierarchy, ln+1, from_t + (to_t - from_t)/2.0, to_t))
    surface_is_moved = true;
   
  // do some coarsening and
  // then update ghost cells through doing refinement if it has finer level
  if(ln < hierarchy->getNumberOfLevels() -1 )
  {
    xfer::CoarsenAlgorithm coarsener(dim);
  

    boost::shared_ptr<xfer::CoarsenSchedule> coarsen_schedule;

    horizon->registerCoarsenActive(coarsener,space_coarsen_op);
      
    coarsen_schedule = coarsener.createSchedule(level, hierarchy->getPatchLevel(ln+1));
    level->getBoxLevel()->getMPI().Barrier();
    coarsen_schedule->coarsenData();

    xfer::RefineAlgorithm post_refiner;

    boost::shared_ptr<xfer::RefineSchedule> refine_schedule;

    

    horizon->registerRKRefinerActive(post_refiner, space_refine_op);
    
    refine_schedule = post_refiner.createSchedule(level, NULL);

    level->getBoxLevel()->getMPI().Barrier();
    refine_schedule->fillData(to_t);
  }

  real_t temp = horizon->maxSurfaceMove(hierarchy, ln, weight_idx);
  if( temp> surface_move_shreshold)
    surface_is_moved = true;
  // copy _a to _p and set _p time to next timestamp

  math::HierarchyCellDataOpsReal<real_t> hcellmath(hierarchy,ln,ln);

  horizon->copyAToP(hcellmath);

  horizon->setLevelTime(level, to_t, to_t);

  level->getBoxLevel()->getMPI().Barrier();

  return surface_is_moved;
}


bool CosmoSim::findHorizon(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  int cnt = 0;
  real_t cur_lambda = 0, delta_lambda = 0;
  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  delta_lambda = AHFinder_dt_frac * pw2(tbox::MathUtilities<double>::Min(
    tbox::MathUtilities<double>::Min(grid_geometry.getDx()[0], grid_geometry.getDx()[1]),
    grid_geometry.getDx()[2]));

  std::cout<<"delta_lambda "<<delta_lambda<<"\n";

  // if we could make sure that the horizon is a shpere
  if(horizon->is_sphere)
  {
    tbox::pout<<"Only checking for sphere horizon!\n";
        return horizon->initSphericalSurface(hierarchy, bssnSim, space_refine_op);
  }
  
  while(cnt < AHFinder_iter_limit)
  {
    initHorizonStep(hierarchy);

    if(!advanceHorizonLevel(hierarchy,
                        0,
                        cur_lambda,
                        cur_lambda + delta_lambda))
    {
      tbox::pout<<"No surface exists or the movement of surface is blow the threshold!\n";
      return false;
    }

    tbox::pout<<cnt<<"\n";
    cnt++;
    cur_lambda += delta_lambda;
  }
  return true;
}

  
/**
 * @brief detect NaNs for field with data_id for in patch
 * 
 * @param patch where we want to fine NaNs
 * @param component id we are looking at 
 */
bool CosmoSim::hasNaNs(
  const boost::shared_ptr<hier::Patch>& patch, idx_t data_id)
{
  boost::shared_ptr<pdat::CellData<double> > d_pdata(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch->getPatchData(data_id)));

  
  boost::shared_ptr<pdat::CellData<double> > w_pdata(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch->getPatchData(weight_idx)));

  
  arr_t d = pdat::ArrayDataAccess::access<DIM, double>(
    d_pdata->getArrayData());
  arr_t w = pdat::ArrayDataAccess::access<DIM, double>(
    w_pdata->getArrayData());

  
  const hier::Box& box = patch->getBox();
  
  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        if(w(i,j,k) > 0 && tbox::MathUtilities< double >::isNaN(d(i,j,k)))
        {
          TBOX_ERROR("NaN detected for variable with id " << data_id <<" "<<i<<" "<<j<<" "<<k<<"\n");
          return 1;
        }
      }
    }
  }
  return 0;
}


/**
 * @brief detect NaNs for all fields in the whole hierarchy
 */
bool CosmoSim::isValid(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy) 
{
  int ln_num = hierarchy->getNumberOfLevels();

  int ln;
  for (ln = 0; ln < ln_num; ln++)
  {
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

      for(int l = 0; l < static_cast<idx_t>(variable_id_list.size()); l++)
      {
        hasNaNs(patch, variable_id_list[l]);
      }
     
    }

  }  // loop over levels
  
  
  return 0;
}


  
} /* namespace cosmo */
