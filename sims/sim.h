#ifndef COSMO_SIM_H
#define COSMO_SIM_H

#include "sim.h"
#include "../ICs/ICs.h"
#include "../cosmo_globals.h"

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/shared_ptr.hpp"
#include <sstream>
#include <iomanip>
#include <cstring>
#include <stdlib.h>

using namespace SAMRAI;

namespace cosmo
{

class CosmoSim
{
protected:
  int step;
  int num_steps;

  std::string simulation_type;
  
  BSSN * bssnSim;


public:
  CosmoSim();
  ~CosmoSim();

  // These functions will be called in main();
  // Each derived class should implement them.
  virtual void init() = 0;
  virtual void runStep() = 0;
  virtual void setICs() = 0;

  void simInit();
  void run();
  void runCommonStepTasks();

  idx_t simNumNaNs();

  boost::shared_ptr<tbox::InputDatabase> d_input_db;
  boost::shared_ptr<tbox::Database> d_cosmo_sim_db;
};

} /* namespace cosmo */

#endif
