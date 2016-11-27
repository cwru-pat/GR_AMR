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

namespace cosmo
{

CosmoSim::CosmoSim()
{
  // Initialize iodata first
  iodata = new IOData(_config["output_dir"]);
  // save a copy of config.txt; print defines
  log_defines(iodata);
  iodata->backupFile(_config.getFileName());

  // fix number of simulation steps
  step = 0;
  num_steps = stoi(_config["steps"]);


  // Store simulation type
  simulation_type = _config["simulation_type"];
}

CosmoSim::~CosmoSim()
{
  std::cout << std::flush;
}

/**
 * @brief      Initialize individual simulation class instances
 */
void CosmoSim::simInit()
{
  // Always use GR fields
  bssnSim = new BSSN();
}

/**
 * @brief      Run the simulation.
 */
void CosmoSim::run()
{
  iodata->log("Running simulation...");

  _timer["loop"].start();
  while(step <= num_steps)
  {
    runStep();
    step++;
  }
  _timer["loop"].stop();

  iodata->log("\nEnding simulation.");
  outputStateInformation();
  iodata->log(_timer.getStateString());
  std::cout << std::flush;
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
  return numNaNs(*bssnSim->fields["DIFFphi_a"]);
}


} /* namespace cosmo */
