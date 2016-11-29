#ifndef COSMO_VACUUM_SIM_H
#define COSMO_VACUUM_SIM_H

#include "sim.h"

namespace cosmo
{

/**
 * derived class based on CosmoSim class (sim.h)
 */
class VacuumSim:
  public CosmoSim,
  public mesh::StandardTagAndInitStrategy,
  public appu::VisDerivedDataStrategy
{
public:
  VacuumSim(){}
  ~VacuumSim(){}

  void init();
  void setICs();
  void initVacuumStep();
  void outputVacuumStep();
  void runVacuumStep();
  void runStep();
  

  virtual bool
   packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& variable_name,
      int depth_id,
      double simulation_time) const;
  
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
};

} /* namespace cosmo */

#endif
