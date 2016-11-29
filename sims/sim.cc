#include "sim.h"
#include "../ICs/ICs.h"
#include "../cosmo_globals.h"

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include "boost/shared_ptr.hpp"
#include <sstream>
#include <iomanip>
#include <cstring>
#include <stdlib.h>

#include <cmath>

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
  std::string simulation_type_in):
  input_db(input_db_in),
  cosmo_sim_db((input_db->getDatabase("CosmoSim")),
  variable_db(hier::VariableDatabase::getDatabase()),
  lstream(l_stream_in),
  dim(dim_in),
  step(0),
  num_steps(cosmo_sim_db->getInteger("steps")),
  simulation_type(simulation_type_in)
{
  t_loop = tbox::TimerManager::getManager()->
    getTimer("loop");
  t_init = tbox::TimerManager::getManager()->
    getTimer("init");
  t_RK_steps = tbox::TimerManager::getManager()->
    getTimer("RK_steps");
}

CosmoSim::~CosmoSim()
{

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
void CosmoSim::run()
{
  iodata->log("Running simulation...");


  t_loop->start();
  while(step <= num_steps)
  {
    runStep();
    step++;
  }
  t_loop->stop();

  tbox::pout<<"\nEnding simulation.";
  outputStateInformation();
}

void CosmoSim::runCommonStepTasks()
{
  // check for NAN every step
  if(simNumNaNs() > 0)
  {
    iodata->log("\nNAN detected!");
    throw 10;
  }

  // progress bar in terminal
  //  io_show_progress(step, num_steps);
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
  return numNaNs(*bssnSim->chi_p);
}


} /* namespace cosmo */
