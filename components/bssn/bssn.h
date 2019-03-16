#ifndef COSMO_BSSN
#define COSMO_BSSN

#include "../../cosmo_includes.h"
#include "BSSNGaugeHandler.h"
#include "../boundaries/sommerfield.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

using namespace SAMRAI;

namespace cosmo
{

/**
 * @brief BSSN Class: evolves BSSN metric fields, computes derived quantities
 */
class BSSN
{
public:
  /* arrays for storing fields */
  BSSN_APPLY_TO_FIELDS(VAR_CREATE)
  BSSN_APPLY_TO_SOURCES(VAR_CREATE)
  BSSN_APPLY_TO_GEN1_EXTRAS(VAR_CREATE)

  std::ostream* lstream;
  std::shared_ptr<tbox::Database>& cosmo_bssn_db;
  const tbox::Dimension& dim;
  real_t KO_damping_coefficient;
  BSSNGaugeHandler * gaugeHandler;
  real_t gd_eta;

  


  
  BSSN(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const tbox::Dimension& dim_in,
    std::shared_ptr<tbox::Database> database_in,
    std::ostream* l_stream_in,
    real_t KO_damping_coefficient_in);

  ~BSSN();

  BSSN_APPLY_TO_FIELDS(RK4_IDX_ALL_CREATE)
  BSSN_APPLY_TO_SOURCES_ARGS(RK4_IDX_CREATE,a)
  BSSN_APPLY_TO_GEN1_EXTRAS_ARGS(RK4_IDX_CREATE,a)
  
  BSSN_APPLY_TO_FIELDS(RK4_PDATA_ALL_CREATE)
  BSSN_APPLY_TO_SOURCES_ARGS(RK4_PDATA_CREATE,a)
  BSSN_APPLY_TO_GEN1_EXTRAS_ARGS(RK4_PDATA_CREATE,a)


  BSSN_APPLY_TO_FIELDS(RK4_MDA_ACCESS_ALL_CREATE)
  BSSN_APPLY_TO_SOURCES_ARGS(RK4_MDA_ACCESS_CREATE,a)
  BSSN_APPLY_TO_GEN1_EXTRAS_ARGS(RK4_MDA_ACCESS_CREATE,a)
  

  void init(const std::shared_ptr<hier::PatchHierarchy>& hierarchy);

  void stepInit(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void RKEvolvePatchBD(const std::shared_ptr<hier::Patch> & patch, real_t dt);

  void set_norm(
    const std::shared_ptr<hier::PatchLevel>& level);
  void set_norm(
    const std::shared_ptr<hier::Patch>& patch, bool need_init_arr);

  void rescale_lapse(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int weight_idx);

  

  void RKEvolvePatch(
    const std::shared_ptr<hier::Patch> & patch, real_t dt);
  void RKEvolvePt(
    idx_t i, idx_t j, idx_t k, BSSNData &bd, const real_t dx[], real_t dt);
  void RKEvolvePtBd(
    idx_t i, idx_t j, idx_t k, BSSNData &bd, const real_t dx[], real_t dt,
    int l_idx, int codim);

  
  void prepareForK1(
    const std::shared_ptr<hier::PatchLevel> & level, real_t to_t);
  void prepareForK2(
    const std::shared_ptr<hier::PatchLevel> & level, real_t to_t);
  void prepareForK3(
    const std::shared_ptr<hier::PatchLevel> & level, real_t to_t);
  void prepareForK4(
    const std::shared_ptr<hier::PatchLevel> & level, real_t to_t);

  void allocField(  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln);
  void allocSrc(  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln);
  void allocGen1(  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln);
  
  void clearSrc(const std::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void clearSrc(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln);
  void clearField(const std::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void clearField(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln);
  void clearGen1(const std::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void clearGen1(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln);

  void addFieldsToList(std::vector<idx_t> &list);
  
  void registerRKRefiner(
    xfer::RefineAlgorithm& refiner,
    std::shared_ptr<hier::RefineOperator> &space_refine_op);
  void registerRKRefinerActive(
    xfer::RefineAlgorithm& refiner,
    std::shared_ptr<hier::RefineOperator> &space_refine_op);

  void registerCoarsenActive(
    xfer::CoarsenAlgorithm& coarsener,
    std::shared_ptr<hier::CoarsenOperator>& coarsen_op);
  void copyAToP(
    math::HierarchyCellDataOpsReal<real_t> & hcellmath);
#if USE_BACKUP_FIELDS
  void copyBToP(
    math::HierarchyCellDataOpsReal<real_t> & hcellmath);
  void copyPToB(
    math::HierarchyCellDataOpsReal<real_t> & hcellmath);
  void copyBToA(
    math::HierarchyCellDataOpsReal<real_t> & hcellmath);
#endif
  void initPData(
    const std::shared_ptr<hier::Patch> & patch);
  void initMDA(
    const std::shared_ptr<hier::Patch> & patch);
  void setLevelTime(
    const std::shared_ptr<hier::PatchLevel> & level,
    double from_t, double to_t);
  void K1FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);
  void K2FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);
  void K3FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);
  void K4FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch, int ln, int max_ln);
  void K4FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);

  void set_bd_values_bd(
    idx_t i, idx_t j, idx_t k, BSSNData *bd, const real_t dx[]);
  void set_bd_values(
    idx_t i, idx_t j, idx_t k, BSSNData *bd, const real_t dx[]);

  void set_local_vals(BSSNData *bd);

  void set_gammai_values(idx_t i, idx_t j, idx_t k, BSSNData *bd);

  void calculate_Acont(BSSNData *bd, const real_t dx[]);

  void calculate_dgamma(BSSNData *bd, const real_t dx[]);

  void calculate_ddgamma(BSSNData *bd, const real_t dx[]);

  void calculate_dalpha_dchi(BSSNData *bd, const real_t dx[]);

  void calculate_dK(BSSNData *bd, const real_t dx[]);
  
#ifdef USE_Z4C
  void calculate_dtheta(BSSNData *bd, const real_t dx[]);
#endif

#ifdef USE_BSSN_SHIFT
  void calculate_dbeta(BSSNData *bd, const real_t dx[]);
#endif

/* #ifdef USE_EXPANSION */
/*   void calculate_dexpN(BSSNData *bd, const real_t dx[]); */
/* #endif */


  
  


    /* Calculate "dependent" quantities (depend on previously calc'd vals) */
  void calculate_conformal_christoffels(BSSNData *bd, const real_t dx[]);

    /* Calculate doubly-"dependent" quantities (depend on previously calc'd vals) */
  void calculateDDchi(BSSNData *bd, const real_t dx[]);
  void calculateRicciTF(BSSNData *bd, const real_t dx[]);
  void calculateDDalphaTF(BSSNData *bd, const real_t dx[]);
  void calculateDZ(BSSNData *bd, const real_t dx[]);
    
  /* Evolution functions */
  real_t ev_DIFFgamma11(BSSNData *bd, const real_t dx[]);
  real_t ev_DIFFgamma12(BSSNData *bd, const real_t dx[]);
  real_t ev_DIFFgamma13(BSSNData *bd, const real_t dx[]);
  real_t ev_DIFFgamma22(BSSNData *bd, const real_t dx[]);
  real_t ev_DIFFgamma23(BSSNData *bd, const real_t dx[]);
  real_t ev_DIFFgamma33(BSSNData *bd, const real_t dx[]);
  real_t ev_A11(BSSNData *bd, const real_t dx[]);
  real_t ev_A12(BSSNData *bd, const real_t dx[]);
  real_t ev_A13(BSSNData *bd, const real_t dx[]);
  real_t ev_A22(BSSNData *bd, const real_t dx[]);
  real_t ev_A23(BSSNData *bd, const real_t dx[]);
  real_t ev_A33(BSSNData *bd, const real_t dx[]);
  real_t ev_DIFFK(BSSNData *bd, const real_t dx[]);
  real_t ev_DIFFchi(BSSNData *bd, const real_t dx[]);
  real_t ev_Gamma1(BSSNData *bd, const real_t dx[]);
  real_t ev_Gamma2(BSSNData *bd, const real_t dx[]);
  real_t ev_Gamma3(BSSNData *bd, const real_t dx[]);

  real_t ev_DIFFalpha(BSSNData *bd, const real_t dx[]);

  real_t ev_theta(BSSNData *bd, const real_t dx[]);
  real_t ev_H(BSSNData *bd, const real_t dx[]);
  real_t ev_a(BSSNData *bd, const real_t dx[]);

#   if USE_BSSN_SHIFT
  real_t ev_beta1(BSSNData *bd, const real_t dx[]);
  real_t ev_beta2(BSSNData *bd, const real_t dx[]);
  real_t ev_beta3(BSSNData *bd, const real_t dx[]);
#   endif

#   if USE_EXPANSION
  real_t ev_expN(BSSNData *bd, const real_t dx[]);
#endif
  
#   if USE_GAMMA_DRIVER
  real_t ev_auxB1(BSSNData *bd, const real_t dx[]);
  real_t ev_auxB2(BSSNData *bd, const real_t dx[]);
  real_t ev_auxB3(BSSNData *bd, const real_t dx[]);
#   endif

#   if USE_PROPER_TIME
  real_t ev_tau(BSSNData *bd, const real_t dx[]);
#   endif

  
  real_t ev_DIFFgamma11_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_DIFFgamma12_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_DIFFgamma13_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_DIFFgamma22_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_DIFFgamma23_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_DIFFgamma33_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_A11_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_A12_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_A13_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_A22_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_A23_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_A33_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_DIFFK_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_DIFFchi_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_Gamma1_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_Gamma2_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_Gamma3_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);

  real_t ev_DIFFalpha_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);

  real_t ev_theta_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_H_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_a_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);

#   if USE_BSSN_SHIFT
  real_t ev_beta1_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_beta2_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_beta3_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
#   endif

#   if USE_EXPANSION
  real_t ev_expN_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
#endif
  
#   if USE_GAMMA_DRIVER
  real_t ev_auxB1_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_auxB2_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
  real_t ev_auxB3_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
#   endif

#   if USE_PROPER_TIME
  real_t ev_tau_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim);
#   endif

  
  void output_max_H_constaint(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    idx_t weight_idx);

  void output_L2_H_constaint(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    idx_t weight_idx, CosmoPatchStrategy * cosmoPS, double exclude_radius);

  void output_L2_H_constaint(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    idx_t weight_idx, CosmoPatchStrategy * cosmoPS);

  /* constraint violation calculations */

  real_t hamiltonianConstraintCalc(BSSNData *bd, const real_t dx[]);
  real_t hamiltonianConstraintScale(BSSNData *bd, const real_t dx[]);

  real_t momentumConstraintCalc(BSSNData *bd, const real_t dx[]);
  real_t momentumConstraintScale(BSSNData *bd, const real_t dx[]);

  void setExtraFieldData();

  void cal_Weyl_scalars(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    int weight_idx);

  // Domain size
  real_t L[DIM];

  bool normalize_Aij, normalize_gammaij;

  real_t Z4c_K1_DAMPING_AMPLITUDE, Z4c_K2_DAMPING_AMPLITUDE;

  real_t chi_lower_bd;

  real_t alpha_lower_bd_for_L2;

  double  domain_lower[3];
  double  domain_upper[3];

  std::vector<double> origin;
  void deBug(  const std::shared_ptr<hier::Patch> & patch);

  double K0;
  double K_avg, rho_P_avg;
};

}

#endif
