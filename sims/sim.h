#ifndef COSMO_SIM_H
#define COSMO_SIM_H

#include "../cosmo_includes.h"
#include "../components/bssn/bssn.h"
#include "../components/boundaries/sommerfield.h"
#include "../cosmo_macros.h"
#include "vacuum.h"

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
  void setGriddingAlgs(
    boost::shared_ptr<mesh::GriddingAlgorithm>& gridding_algorithm_in);

  idx_t simNumNaNs();

  boost::shared_ptr<tbox::InputDatabase>& input_db;
  boost::shared_ptr<tbox::Database>& cosmo_sim_db;
  hier::VariableDatabase* variable_db;
  boost::shared_ptr<mesh::GriddingAlgorithm>& gridding_algorithm;
  std::ostream* lstream;

  const tbox::Dimension& dim;

  std::string simulation_type;

  idx_t do_plot;
  real_t dt_frac;

  std::string vis_filename;

  real_t cur_t;
  static boost::shared_ptr<tbox::Timer> t_loop;
  static boost::shared_ptr<tbox::Timer> t_init;
  static boost::shared_ptr<tbox::Timer> t_RK_steps;
  // patch strategy that managers boundary and refine
  // && coarsen strategy
  CosmoPatchStrategy * cosmoPS;

  CosmoIO *cosmo_io;

  idx_t weight_idx;
  boost::shared_ptr<pdat::CellVariable<real_t> > weight;
  idx_t regridding_interval;
};

} /* namespace cosmo */

#endif
