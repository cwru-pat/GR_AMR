#ifndef COSMO_SCALAR
#define COSMO_SCALAR

#include "../../cosmo_includes.h"
#include "scalar_macros.h"
#include "../boundaries/sommerfield.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "scalarPotentialHandler.h"
#include "scalar_data.h"
#include "../../components/bssn/bssn.h"

using namespace SAMRAI;


namespace cosmo
{
 
class Scalar
{
public:

  std::ostream* lstream;

  std::shared_ptr<tbox::Database>& cosmo_scalar_db;

  const tbox::Dimension& dim;

  Scalar(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const tbox::Dimension& dim_in,
    std::shared_ptr<tbox::Database> database_in,
    std::ostream* l_stream_in,
    real_t KO_damping_coefficient_in);

  ~Scalar();

  SCALAR_APPLY_TO_FIELDS(VAR_CREATE)
  
  SCALAR_APPLY_TO_FIELDS(RK4_IDX_ALL_CREATE)
  SCALAR_APPLY_TO_FIELDS(RK4_PDATA_ALL_CREATE)
  SCALAR_APPLY_TO_FIELDS(RK4_MDA_ACCESS_ALL_CREATE)

  real_t KO_damping_coefficient;

  scalarPotentialHandler * potentialHandler;
  
  void getScalarData(
    idx_t i, idx_t j, idx_t k, BSSNData *bd, ScalarData *sd, const real_t dx[]);  
  void getScalarDataBd(
    idx_t i, idx_t j, idx_t k, BSSNData *bd, ScalarData * sd, const real_t dx[]);

  
  void alloc(const std::shared_ptr<hier::PatchHierarchy> &hierarchy, idx_t ln);
  
  void addFieldsToList(std::vector<idx_t> &list);

  void registerRKRefinerActive(
    xfer::RefineAlgorithm& refiner,
    std::shared_ptr<hier::RefineOperator> &space_refine_op);

  
  void registerRKRefiner(
    xfer::RefineAlgorithm& refiner,
    std::shared_ptr<hier::RefineOperator> &space_refine_op);
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
  void addBSSNSrc(
    BSSN * bssn, const std::shared_ptr<hier::Patch> & patch, bool need_init_arr);
  void addBSSNSrc(
    BSSN * bssn, const std::shared_ptr<hier::PatchLevel> & level);
  void addBSSNSrc(
    BSSN * bssn, const std::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void clear(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln);
  void registerCoarsenActive(
    xfer::CoarsenAlgorithm& coarsener,
    std::shared_ptr<hier::CoarsenOperator>& coarsen_op);

  
  void K1FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);
  void K2FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);
  void K3FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);
  void K4FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);
  void init(const std::shared_ptr<hier::PatchHierarchy>& hierarchy);

  void stepInit(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy);
  
  void RKEvolvePatchBD(const std::shared_ptr<hier::Patch> & patch, real_t dt);

  void set_norm(
    const std::shared_ptr<hier::PatchLevel>& level);

  

  void RKEvolvePatch(
    const std::shared_ptr<hier::Patch> & patch, real_t dt);

  void RKEvolvePt(
    idx_t i, idx_t j, idx_t k, BSSNData &bd, ScalarData &sd, const real_t dx[], real_t dt);

  void RKEvolvePtBd(
    idx_t i, idx_t j, idx_t k, BSSNData &bd, ScalarData &sd,
    const real_t dx[], real_t dt, int l_idx, int codim);

  
  void prepareForK1(
    const std::shared_ptr<hier::PatchLevel> & level, real_t to_t);
  void prepareForK2(
    const std::shared_ptr<hier::PatchLevel> & level, real_t to_t);
  void prepareForK3(
    const std::shared_ptr<hier::PatchLevel> & level, real_t to_t);
  void prepareForK4(
    const std::shared_ptr<hier::PatchLevel> & level, real_t to_t);

  real_t ev_phi(BSSNData *bd, ScalarData *sd, const real_t dx[]);
  real_t ev_Pi(BSSNData *bd, ScalarData *sd, const real_t dx[]);  
  real_t ev_psi1(BSSNData *bd, ScalarData *sd, const real_t dx[]);
  real_t ev_psi2(BSSNData *bd, ScalarData *sd, const real_t dx[]);
  real_t ev_psi3(BSSNData *bd, ScalarData *sd, const real_t dx[]);

  real_t ev_phi_bd(BSSNData *bd, ScalarData *sd, const real_t dx[], int l_idx, int codim);
  real_t ev_Pi_bd(BSSNData *bd, ScalarData *sd, const real_t dx[], int l_idx, int codim);  
  real_t ev_psi1_bd(BSSNData *bd, ScalarData *sd, const real_t dx[], int l_idx, int codim);
  real_t ev_psi2_bd(BSSNData *bd, ScalarData *sd, const real_t dx[], int l_idx, int codim);
  real_t ev_psi3_bd(BSSNData *bd, ScalarData *sd, const real_t dx[], int l_idx, int codim);

  
  
  


  

};

} // namespace cosmo


#endif
