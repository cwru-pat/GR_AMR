#ifndef COSMO_SIM_H
#define COSMO_SIM_H

#include "../cosmo_includes.h"
#include "../components/bssn/bssn.h"
#include "../components/IO/io.h"
#include "../cosmo_ps.h"
#include "../cosmo_macros.h"
#include "../cosmo_types.h"

using namespace SAMRAI;

namespace cosmo
{

class CosmoSim:public mesh::StandardTagAndInitStrategy
{  
public:

  BSSN * bssnSim;
  CosmoSim(
    const tbox::Dimension& dim_in,
    boost::shared_ptr<tbox::InputDatabase>& input_db_in,
    std::ostream* l_stream_in,
    std::string simulation_type_in,
    std::string vis_filename_in);

  virtual ~CosmoSim(
    void);


  // These functions will be called in main();
  // Each derived class should implement them.
  virtual void init() = 0;
  virtual void runStep(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy) =0;
  virtual void setICs() = 0;

    virtual void
    initializeLevelData(
      /*! Hierarchy to initialize */
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      /*! Level to initialize */
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      /*! Whether level is being introduced for the first time */
      const bool initial_time,
      /*! Level to copy data from */
      const boost::shared_ptr<hier::PatchLevel>& old_level =
      boost::shared_ptr<hier::PatchLevel>(),
      /*! Whether data on new patch needs to be allocated */
      const bool allocate_data = true) = 0;

  virtual void
    resetHierarchyConfiguration(
      /*! New hierarchy */
      const boost::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
      /*! Coarsest level */ int coarsest_level,
      /*! Finest level */ int finest_level) = 0;

  virtual void
    applyGradientDetector(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double error_data_time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation) = 0;

  
  void setRefineCoarsenOps(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void run(  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void runCommonStepTasks(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void setGriddingAlgs(
    boost::shared_ptr<mesh::GriddingAlgorithm>& gridding_algorithm_in);

  bool isValid(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);
  bool hasNaNs(
    const boost::shared_ptr<hier::Patch>& patch, idx_t data_id);

  
  boost::shared_ptr<tbox::InputDatabase>& input_db;
  boost::shared_ptr<tbox::Database> cosmo_sim_db;
  //hier::VariableDatabase* variable_db;
  boost::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm;
  std::ostream* lstream;

  const tbox::Dimension dim;

  int step;
  int num_steps;

  std::string simulation_type;

  bool do_plot;
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

  boost::shared_ptr<pdat::CellVariable<real_t> > weight;
  idx_t weight_idx;

  idx_t regridding_interval;
  real_t KO_damping_coefficient;
  real_t adaption_threshold;
  boost::shared_ptr<hier::RefineOperator> space_refine_op;
  boost::shared_ptr<hier::CoarsenOperator> space_coarsen_op;

  std::vector<int> variable_id_list;
  
};

} /* namespace cosmo */

#endif
