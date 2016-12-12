#ifndef COSMO_BSSN
#define COSMO_BSSN

#include "../../cosmo_includes.h"
#include "bssn.h"
#include "BSSNGaugeHandler.h"
#include "../boundaries/sommerfield.h"

namespace cosmo
{

/**
 * @brief BSSN Class: evolves BSSN metric fields, computes derived quantities
 */
class BSSN
{
  /* arrays for storing fields */
  BSSN_APPLY_TO_FIELDS(VAR_CREATE);
  BSSN_APPLY_TO_SOURCES(VAR_CREATE)
  BSSN_APPLY_TO_GEN1_EXTRAS(VAR_CREATE)

  BSSNGaugeHandler * gaugeHandler;

  real_t KO_damping_coefficient;
  real_t gd_eta;

  boost::shared_ptr<tbox::Database> d_bssn_db;
  
public:

  BSSN(
  const tbox::Dimension& dim_in,
  tbox::Database& database_in,
  std::ostream* l_stream_in = 0
  xfer::RefinePatchStrategy* PS_in);

  ~BSSN();

  BSSN_APPLY_TO_FIELDS(RK4_IDX_ALL_CREATE);
  BSSN_APPLY_TO_SOURCES_ARGS(RK4_IDX_CREATE,a);
  BSSN_APPLY_TO_GEN1_EXTRAS_ARGS(RK4_IDX_CREATE,a);
  
  BSSN_APPLY_TO_FIELDS(RK4_PDATA_ALL_CREATE);
  BSSN_APPLY_TO_SOURCES(RK4_PDATA_CREATE);
  BSSN_APPLY_TO_GEN1_EXTRAS(RK4_PDATA_CREATE);


  BSSN_APPLY_TO_FIELDS(RK4_MDA_ACCESS_ALL_CREATE);
  BSSN_APPLY_TO_SOURCES_ARGS(RK4_MDA_ACCESS_CREATE,a);
  BSSN_APPLY_TO_GEN1_EXTRAS_ARGS(RK4_MDA_ACCESS_CREATE,a);  
  
  boost::shared_ptr<hier::RefineOperator> space_refine_op;
  boost::shared_ptr<hier::CoarsenOperator> space_coarsen_op;

  void init();

  void RKEvolvePatchBD(const boost::shared_ptr<hier::Patch> & patch, real_t dt);

#if USE_CCZ4
  void initZ(
    const boost::shared_ptr<hier::PatchLevel> & level);
#endif

  void RKEvolvePatch(
    const boost::shared_ptr<hier::Patch> & patch, real_t dt);

  void prepairForK1(
    const boost::shared_ptr<hier::PatchLevel> & level, real_t to_t);
  void prepairForK2(
    const boost::shared_ptr<hier::PatchLevel> & level, real_t to_t);
  void prepairForK3(
    const boost::shared_ptr<hier::PatchLevel> & level, real_t to_t);
  void prepairForK4(
    const boost::shared_ptr<hier::PatchLevel> & level, real_t to_t);

  void registerRKRefiner(
    xfer::RefineAlgorithm& refiner,
    boost::shared_ptr<hier::RefineOperator> &space_refine_op);
  void registerSameLevelRefinerActive(
    xfer::RefineAlgorithm& refiner,
    boost::shared_ptr<hier::RefineOperator> &space_refine_op);
  void registerCoarsenActive(
    xfer::CoarsenAlgorithm& coarsener,
    boost::shared_ptr<hier::CoarsenOperator>& coarsen_op);
  void swapPF(
    math::HierarchyCellDataOpsReal<real_t> & hcellmath);
  void copyPToA(
    math::HierarchyCellDataOpsReal<real_t> & hcellmath);
  void initPData(
    const boost::shared_ptr<hier::Patch> & patch);
  void initMDA(
    const boost::shared_ptr<hier::Patch> & patch);
  void setLevelTime(
    const boost::shared_ptr<hier::PatchLevel> & level,
    double from_t, double to_t);
  void K1FinalizePatch(
    const boost::shared_ptr<hier::Patch> & patch);
  void K2FinalizePatch(
    const boost::shared_ptr<hier::Patch> & patch);
  void K3FinalizePatch(
    const boost::shared_ptr<hier::Patch> & patch);
  void K4FinalizePatch(
    const boost::shared_ptr<hier::Patch> & patch);

  void BSSN::set_bd_values_bd(
    idx_t i, idx_t j, idx_t k, BSSNData *bd, real_t dx[]);
  void BSSN::set_bd_values_for_extra_fields(
    idx_t i, idx_t j, idx_t k, BSSNData *bd, real_t dx[]);
  void BSSN::set_bd_values(
    idx_t i, idx_t j, idx_t k, BSSNData *bd, real_t dx[]);

  void BSSN::set_local_vals(BSSNData *bd);

  void BSSN::set_gammai_values(idx_t i, idx_t j, idx_t k, BSSNData *bd);

  void BSSN::calculate_Acont(BSSNData *bd, real_t dx[]);

  void BSSN::calculate_dgamma(BSSNData *bd, real_t dx[]);

  void BSSN::calculate_ddgamma(BSSNData *bd, real_t dx[]);

  void BSSN::calculate_dalpha_dphi(BSSNData *bd, real_t dx[]);

  void BSSN::calculate_dK(BSSNData *bd, real_t dx[]);
  
#ifdef USE_CCZ4
  void BSSN::calculate_dtheta(BSSNData *bd, real_t dx[]);
#endif

#ifdef USE_BSSN_SHIFT
  void BSSN::calculate_dbeta(BSSNData *bd, real_t dx[]);
#endif

#ifdef USE_EXPANSION
  void BSSN::calculate_dexpN(BSSNData *bd, real_t dx[]);
#endif


  
  


    /* Calculate "dependent" quantities (depend on previously calc'd vals) */
  void calculate_conformal_christoffels(BSSNData *bd, real_t dx[]);

    /* Calculate doubly-"dependent" quantities (depend on previously calc'd vals) */
  void calculateDDphi(BSSNData *bd, real_t dx[]);
  void calculateRicciTF(BSSNData *bd, real_t dx[]);
  void calculateDDalphaTF(BSSNData *bd, real_t dx[]);
  void BSSN::calculateDZ(BSSNData *bd, real_t dx[]);
    
  /* Evolution functions */
  real_t ev_DIFFgamma11(BSSNData *bd, real_t dx[]);
  real_t ev_DIFFgamma12(BSSNData *bd, real_t dx[]);
  real_t ev_DIFFgamma13(BSSNData *bd, real_t dx[]);
  real_t ev_DIFFgamma22(BSSNData *bd, real_t dx[]);
  real_t ev_DIFFgamma23(BSSNData *bd, real_t dx[]);
  real_t ev_DIFFgamma33(BSSNData *bd, real_t dx[]);
  real_t ev_A11(BSSNData *bd, real_t dx[]);
  real_t ev_A12(BSSNData *bd, real_t dx[]);
  real_t ev_A13(BSSNData *bd, real_t dx[]);
  real_t ev_A22(BSSNData *bd, real_t dx[]);
  real_t ev_A23(BSSNData *bd, real_t dx[]);
  real_t ev_A33(BSSNData *bd, real_t dx[]);
  real_t ev_DIFFK(BSSNData *bd, real_t dx[]);
  real_t ev_DIFFphi(BSSNData *bd, real_t dx[]);
  real_t ev_Gamma1(BSSNData *bd, real_t dx[]);
  real_t ev_Gamma2(BSSNData *bd, real_t dx[]);
  real_t ev_Gamma3(BSSNData *bd, real_t dx[]);

  real_t ev_DIFFalpha(BSSNData *bd, real_t dx[]);

  real_t ev_theta(BSSNData *bd, real_t dx[]);

#   if USE_BSSN_SHIFT
  real_t ev_beta1(BSSNData *bd, real_t dx[]);
  real_t ev_beta2(BSSNData *bd, real_t dx[]);
  real_t ev_beta3(BSSNData *bd, real_t dx[]);
#   endif

#   if USE_EXPANSION
  real_t ev_expN(BSSNData *bd, real_t dx[]);
#endif
  
#   if USE_GAMMA_DRIVER
  real_t ev_auxB1(BSSNData *bd, real_t dx[]);
  real_t ev_auxB2(BSSNData *bd, real_t dx[]);
  real_t ev_auxB3(BSSNData *bd, real_t dx[]);
#   endif

  void BSSN::output_max_H_constaint(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy
    idx_t weight_idx);

  /* constraint violation calculations */

  real_t hamiltonianConstraintCalc(BSSNData *bd, real_t dx[]);
  real_t hamiltonianConstraintScale(BSSNData *bd, real_t dx[]);



};

}

#endif
