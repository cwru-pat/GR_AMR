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
  use_anguler_momentum_finder(cosmo_sim_db->getBoolWithDefault("use_anguler_momentum_finder", false)),
  gradient_indicator(cosmo_sim_db->getStringWithDefault("gradient_indicator", "DIFFchi")),
  regridding_step_lower_bound(cosmo_sim_db->getIntegerWithDefault("regridding_step_lower_bound", 0)),
  regridding_step_upper_bound(cosmo_sim_db->getIntegerWithDefault("regridding_step_upper_bound", 100000000)),
  stop_after_found_horizon(cosmo_sim_db->getBoolWithDefault("stop_after_found_horizon",false))
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

  horizon = new AHFinderDirect::Horizon(
    hierarchy, bssnSim, dim,input_db->getDatabase("AHFD"), vis_filename.c_str(),weight_idx);

  horizon->AHFinderDirect_setup();

  // scractch component for refinement, currently not in use
  refine_scratch_idx = variable_db->registerVariableAndContext( 
      refine_scratch, 
      context_active,
      hier::IntVector(dim, STENCIL_ORDER));

  AHFinder_interval = cosmo_sim_db->getIntegerWithDefault("AHFinder_interval" , 0);

  if(cosmo_sim_db->keyExists("save_steps"))
    save_steps = cosmo_sim_db->getIntegerVector("save_steps");
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
bool CosmoSim::runCommonStepTasks(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
    // detecting NaNs
  isValid(hierarchy);

  if(step > starting_step &&
    ((step %save_interval == 0) ||
      (std::find(save_steps.begin(), save_steps.end(), step) != save_steps.end()) ))
  {
    std::string restart_file_name = simulation_type + comments +".restart";
    tbox::RestartManager::getManager()->writeRestartFile(restart_file_name, step);
  }

  // since all neccecery levels were built when setting IC
  // no need to regrid again at zero step
  if(step > starting_step && step >= regridding_step_lower_bound
     && step <= regridding_step_upper_bound
     && (step % regridding_interval == 0))
  {
    std::vector<int> tag_buffer(hierarchy->getMaxNumberOfLevels());
    for (idx_t ln = 0; ln < static_cast<int>(tag_buffer.size()); ++ln) {
      tag_buffer[ln] = 1;
    }
    gridding_algorithm->regridAllFinerLevels(
      0,
      tag_buffer,
      step,
      cur_t);
    tbox::plog << "Newly adapted hierarchy\n";
    hierarchy->recursivePrint(tbox::plog, "    ", 1);
  }

  bool found_horizon = false;
  
  if(use_AHFinder)
  {
    horizon->AHFinderDirect_find_horizons(step, cur_t);
    //    found_horizon = horizon->AHFinderDirect_horizon_was_found(1);
  }
  // if(found_horizon && use_anguler_momentum_finder)
  // {
  //   horizon->findKilling(hierarchy, bssnSim);
    
  // }
  return found_horizon;
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

#pragma omp parallel for collapse(2)
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        if(w(i,j,k) > 0 && tbox::MathUtilities< double >::isNaN(d(i,j,k)))
        {
          TBOX_ERROR("NaN detected for variable with id " << data_id <<" "<<i<<" "<<j<<" "<<k<<"\n");
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
