#ifndef COSMO_VACUUM_SIM_H
#define COSMO_VACUUM_SIM_H

#include "sim.h"
#include "../cosmo_includes.h"
#include "../components/boundaries/sommerfield.h"
#include "../components/boundaries/periodic.h"
#include "../components/bssn/bssn.h"
#include "../components/bssn/bssn_ic.h"
#include "vacuum_macros.h"


using namespace SAMRAI;

namespace cosmo
{

/**
 * derived class based on CosmoSim class (sim.h)
 */
class VacuumSim:
  public CosmoSim  
{
public:
  VacuumSim(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const tbox::Dimension& dim_in,
    boost::shared_ptr<tbox::InputDatabase>& input_db_in,
    std::ostream* l_stream_in,
    std::string simulation_type_in,
    std::string vis_filename_in);

  virtual ~VacuumSim(
    void);


  void init();
  void setICs(const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void initVacuumStep(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void outputVacuumStep(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void runVacuumStep(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double from_t, double to_t);
  virtual void runStep(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);
  double getDt(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);
  

  
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
      const bool allocate_data = true);

  virtual void
    resetHierarchyConfiguration(
      /*! New hierarchy */
      const boost::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
      /*! Coarsest level */ int coarsest_level,
      /*! Finest level */ int finest_level);

  virtual void
    applyGradientDetector(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double error_data_time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation);
 private:
  boost::shared_ptr<tbox::Database> cosmo_vacuum_db;

  void initLevel(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln);
  void computeVectorWeights(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void addBSSNExtras(
    const boost::shared_ptr<hier::PatchLevel> & level);
  void addBSSNExtras(
    const boost::shared_ptr<hier::Patch> & patch);
  void RKEvolveLevel(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    idx_t ln,
    double from_t,
    double to_t);
  void advanceLevel(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    int ln,
    double from_t,
    double to_t);

  
  
};

} /* namespace cosmo */

#endif
