
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/solv/CellPoissonFACOps.h"
#include "SAMRAI/tbox/Dimension.h"


#include <string>

#include "SAMRAI/tbox/Database.h"

/*
 * SAMRAI classes
 */
#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/geom/CartesianCellDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianSideDoubleWeightedAverage.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/solv/CartesianRobinBcHelper.h"
#include "SAMRAI/solv/FACPreconditioner.h"
#include "SAMRAI/solv/GhostCellRobinBcCoefs.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include "boost/shared_ptr.hpp"

using namespace SAMRAI;

class WaveSolv:
   public mesh::StandardTagAndInitStrategy,
   public appu::VisDerivedDataStrategy
{

 public:
  WaveSolv(const tbox::Dimension& dim,
           tbox::Database& database,
           std::ostream* log_stream);

  double _laplacian(  MDA_Access<double, 3, MDA_OrderColMajor<3> > &var, int i, int j, int k, const double dx[]);
  double _der_norm(boost::shared_ptr<pdat::CellData<double>>& var, int i, int j, int k, const double dx[]);
  double _maxError(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);

  void initCoarsest(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);  

  void computeVectorWeights(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);



  void evolve_patch(
   double* phi_p,
   double* pi_p,
   double* phi_c,
   double* pi_c,
   const int phigi,
   const int phigj,
   const int phigk,
   const int pigi,
   const int pigj,
   const int pigk,
   const int ifirst,
   const int ilast,
   const int jfirst,
   const int jlast,
   const int kfirst,
   const int klast,
   const double dx[],
   const double dt);

  
  void advanceLevel(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  int ln,
  double from_t,
  double to_t);

  void advanceHierarchy(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double from_t, double to_t);

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

  std::ostream* d_lstream;
  tbox::Database *d_database;
  bool d_barrier_and_time;
  static boost::shared_ptr<tbox::Timer> t_advance_hier;

  double d_omega;
  double nx, ny, nz;
  /*#ifdef HAVE_HDF5
  boost::shared_ptr<appu::VisItDataWriter> d_visit_writer;
  #endif*/

#ifdef HAVE_HDF5
   /*!
    * @brief Tell a plotter which data to write for this class.
    */
   int
   registerVariablesWithPlotter(
      appu::VisItDataWriter& visit_writer);
#endif


  
 public:
  const tbox::Dimension d_dim;

  boost::shared_ptr<hier::PatchHierarchy> d_hierarchy;

  boost::shared_ptr<hier::VariableContext> d_context_current;

  boost::shared_ptr<hier::VariableContext> d_context_scratch;
  boost::shared_ptr<hier::VariableContext> d_context_previous;

  boost::shared_ptr<pdat::CellVariable<double> > d_phi;
  boost::shared_ptr<pdat::CellVariable<double> > d_pi;
  boost::shared_ptr<pdat::CellVariable<double> > d_wt;

  int d_phi_current, d_phi_scratch, d_phi_previous;

  int d_pi_current, d_pi_scratch, d_pi_previous;

  int d_weight;
  
  double d_adaption_threshold;

#ifdef HAVE_HDF5
   boost::shared_ptr<appu::VisItDataWriter> d_visit_writer;
#endif
   int d_finest_dbg_plot_ln;
   //@}  
  

};
