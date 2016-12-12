#include "sim.h"
#include "../cosmo_includes.h"
#include "../cosmo_macros.h"
#include "../components/bssn/bssn.h"

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/MathUtilities.h"


using namespace SAMRAI;

boost::shared_ptr<tbox::Timer> CosmoSim::t_loop;

boost::shared_ptr<tbox::Timer> VacuumSim::t_init;
boost::shared_ptr<tbox::Timer> VacuumSim::t_RK_steps;

namespace cosmo
{

CosmoSim::CosmoSim(
  const tbox::Dimension& dim_in,
  boost::shared_ptr<tbox::InputDatabase>& input_db_in,
  std::ostream* l_stream_in = 0,
  std::string simulation_type_in,
  std::string vis_filename_in):
  input_db(input_db_in),
  cosmo_sim_db((input_db->getDatabase("CosmoSim")),
  variable_db(hier::VariableDatabase::getDatabase()),
  gridding_algorithm(gridding_algorithm_in),
  lstream(l_stream_in),
  dim(dim_in),
  step(0),
  num_steps(cosmo_sim_db->getInteger("steps")),
  simulation_type(simulation_type_in),
  do_plot(cosmo_sim_db->getBoolWithDefault("do_plot", false)),
  dt_frac(cosmo_sim_db->getDoubleWithDefault("dt_frac", 0.1)),
  vis_filename(vis_filename_in),
  cur_t(0),
  weight(new pdat::CellVariable<real_t>(dim, "weight", 1)),
  weight_idx(0),
  regridding_interval(cosmo_sim_db->getInteger("regridding_interval"))
{
  t_loop = tbox::TimerManager::getManager()->
    getTimer("loop");
  t_init = tbox::TimerManager::getManager()->
    getTimer("init");
  t_RK_steps = tbox::TimerManager::getManager()->
    getTimer("RK_steps");
  simInit();
  cosmo_io = new CosmoIO(dim, input_db->getDataBase("IO"), lstream);

  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  
  boost::shared_ptr<hier::VariableContext> context_active(
    variable_db->getContext("ACTIVE"));

  weight_idx = variable_db->registerVariableAndContext( 
      weight, 
      context,
      hier::IntVector(dim, STENCIL_ORDER));
}

CosmoSim::~CosmoSim()
{

}

void CosmoSim::setGriddingAlgs(
  boost::shared_ptr<mesh::GriddingAlgorithm>& gridding_algorithm_in)
{
  gridding_algorithm = gridding_algorithm_in;
}
/**
 * @brief      Initialize individual simulation class instances
 */
void CosmoSim::simInit()
{
  // Always use GR fields
  bssnSim = new BSSN(dim,input_db->getDatabase("BSSN"), lstream);
}

/**
 * @brief      Run the simulation.
 */
void CosmoSim::run(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  tbox::plog("Running simulation...");

  t_loop->start();
  while(step <= num_steps)
  {
    runStep(hierarchy);
    step++;
  }
  t_loop->stop();

  tbox::pout<<"\nEnding simulation.";
  outputStateInformation();
}

               
void CosmoSim::runCommonStepTasks(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  if(step % regridding_interval == 0)
  {
    std::vector<int> tag_buffer(hierarchy->getMaxNumberOfLevels());
    for (ln = 0; ln < static_cast<int>(tag_buffer.size()); ++ln) {
      tag_buffer[ln] = 1;
    }
    gridding_algorithm->regridAllFinerLevels(
      0,
      tag_buffer,
      0,
      0.0);
    tbox::plog << "Newly adapted hierarchy\n";
    hierarchy->recursivePrint(tbox::plog, "    ", 1);
  }
}

  //TODO
void CosmoSim::prepBSSNOutput()
{
}

  //TODO
void CosmoSim::outputStateInformation()
{
  return;
}

  //TODO
int CosmoSim::simNumNaNs() 
{
  // check for NAN in a field
  //return numNaNs(*bssnSim->DIFFchi_p);
}


} /* namespace cosmo */
