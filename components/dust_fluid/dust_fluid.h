#ifndef COSMO_DUST_FLUID
#define COSMO_DUST_FLUID

#include "../../cosmo_includes.h"
#include "../bssn/bssn.h"
#include "../bssn/bssn_ic.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "dust_fluid_macros.h"
#include "dust_fluid_data.h"

namespace cosmo
{

/** Static matter class **/
class DustFluid
{
  /* Fluid field */
  // just a density variable
public:

  std::ostream* lstream;
  std::shared_ptr<tbox::Database>& cosmo_dust_fluid_db;
  const tbox::Dimension& dim;
  
  DustFluid(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const tbox::Dimension& dim_in,
    std::shared_ptr<tbox::Database> database_in,
    std::ostream* l_stream_in);
  
  ~DustFluid();
  
  DUST_FLUID_APPLY_TO_FIELDS(VAR_CREATE)
  DUST_FLUID_APPLY_TO_DERIVED_FIELDS(VAR_CREATE)
  
  DUST_FLUID_APPLY_TO_FIELDS(RK4_IDX_ALL_CREATE)
  DUST_FLUID_APPLY_TO_FIELDS(RK4_PDATA_ALL_CREATE)
  DUST_FLUID_APPLY_TO_FIELDS(RK4_MDA_ACCESS_ALL_CREATE)

  DUST_FLUID_APPLY_TO_DERIVED_FIELDS_ARGS(RK4_IDX_CREATE,a)
  DUST_FLUID_APPLY_TO_DERIVED_FIELDS_ARGS(RK4_PDATA_CREATE,a)
  DUST_FLUID_APPLY_TO_DERIVED_FIELDS_ARGS(RK4_MDA_ACCESS_CREATE,a)
  
  void init(const std::shared_ptr<hier::PatchHierarchy>& hierarchy);

  void alloc(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln);
  
  void clear(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln);
  void clearDerivedFields(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln);

  void stepInit(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy);

  void initPData(
    const std::shared_ptr<hier::Patch> & patch);
  void initMDA(
    const std::shared_ptr<hier::Patch> & patch);

  
  void addBSSNSrc(
    BSSN *bssn,   const std::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void addBSSNSrc(
    BSSN *bssn, const std::shared_ptr<hier::PatchLevel> & level);
  void addBSSNSrc(
    BSSN *bssn, const std::shared_ptr<hier::Patch> & patch, bool need_init_arr);


  void registerRKRefinerActive(
    xfer::RefineAlgorithm& refiner,
    std::shared_ptr<hier::RefineOperator> &space_refine_op);

  void registerRKRefiner(
    xfer::RefineAlgorithm& refiner,
    std::shared_ptr<hier::RefineOperator> &space_refine_op);

  void registerCoarsenActive(
    xfer::CoarsenAlgorithm& coarsener,
    std::shared_ptr<hier::CoarsenOperator>& coarsen_op);

  void copyAToP(
    math::HierarchyCellDataOpsReal<real_t> & hcellmath);
  
  void addFieldsToList(std::vector<idx_t> &list);

  void addDerivedFields(
    BSSN *bssn,   const std::shared_ptr<hier::PatchHierarchy>& hierarchy);  
  void addDerivedFields(
    BSSN *bssn, const std::shared_ptr<hier::PatchLevel> & level);
  void addDerivedFields(
    BSSN *bssn, const std::shared_ptr<hier::Patch> & patch );

  void dust_fluid_ic_set_fluid_for_BHL(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    idx_t ln, 
    std::shared_ptr<tbox::Database> cosmo_dust_fluid_db);

  void prepareForK1(
  const std::shared_ptr<hier::PatchLevel> & level,
  real_t to_t);
  void prepareForK2(
  const std::shared_ptr<hier::PatchLevel> & level,
  real_t to_t);
  void prepareForK3(
  const std::shared_ptr<hier::PatchLevel> & level,
  real_t to_t);
  void prepareForK4(
  const std::shared_ptr<hier::PatchLevel> & level,
  real_t to_t);

  void K1FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);
  void K2FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);
  void K3FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);
  void K4FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);

  real_t ev_DF_D(BSSNData *bd, DustFluidData *dd, const real_t dx[]);
  real_t ev_DF_S1(BSSNData *bd, DustFluidData *dd, const real_t dx[]);
  real_t ev_DF_S2(BSSNData *bd, DustFluidData *dd, const real_t dx[]);
  real_t ev_DF_S3(BSSNData *bd, DustFluidData *dd, const real_t dx[]);
  real_t ev_DF_E(BSSNData *bd, DustFluidData *dd, const real_t dx[]);

  real_t ev_DF_D_bd(BSSNData *bd, DustFluidData *sd, const real_t dx[], int l_idx, int codim);
  real_t ev_DF_S1_bd(BSSNData *bd, DustFluidData *sd, const real_t dx[], int l_idx, int codim);
  real_t ev_DF_S2_bd(BSSNData *bd, DustFluidData *sd, const real_t dx[], int l_idx, int codim);
  real_t ev_DF_S3_bd(BSSNData *bd, DustFluidData *sd, const real_t dx[], int l_idx, int codim);
  real_t ev_DF_E_bd(BSSNData *bd, DustFluidData *sd, const real_t dx[], int l_idx, int codim);


  
  void RKEvolvePatchBD(
  const std::shared_ptr<hier::Patch> & patch,
  real_t dt);

  void RKEvolvePatch(
    const std::shared_ptr<hier::Patch> & patch, BSSN *bssn, real_t dt);
  void RKEvolvePt(
    idx_t i, idx_t j, idx_t k, BSSNData &bd, DustFluidData & dd, const real_t dx[], real_t dt);
  void RKEvolvePtBd(
    idx_t i, idx_t j, idx_t k, BSSNData &bd, DustFluidData &dd,
    const real_t dx[], real_t dt, int l_idx, int codim);

  
  void setLevelTime(
    const std::shared_ptr<hier::PatchLevel> & level,
    double from_t, double to_t);
  void getDustFluidData(
    idx_t i, idx_t j, idx_t k, BSSNData *bd, DustFluidData * dd, const real_t dx[]);
  void getDustFluidDataBd(
    idx_t i, idx_t j, idx_t k, BSSNData *bd, DustFluidData * dd, const real_t dx[]);

  void printWConstraint(
    BSSN *bssn,   const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t weight_idx);


  
  // if true, will not change T_{\mu\nu}
  bool is_test_fluid;
  
};

}

#endif
