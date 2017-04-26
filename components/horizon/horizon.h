#ifndef COSMO_HORIZON
#define COSMO_HORIZON

#include "../../cosmo_includes.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "horizon_data.h"
#include "../components/bssn/bssn.h"

using namespace SAMRAI;


#define d2d1F d1d2F
#define d3d1F d1d3F
#define d3d2F d2d3F


#define HORIZON_CALCULATE_DS_TERM1(J, L, I, K) \
  (bd->gammai##I##J * hd->d##J##F * bd->gammai##K##L * hd->d##L##F)

#define HORIZON_CALCULATE_DS_TERM2(I, K) \
  (hd->d##I##d##K##F - bd->G##1##I##K * hd->d1F   \
   - bd->G##2##I##K * hd->d2F                     \
   - bd->G##3##I##K * hd->d3F                     \
   + (bd->d##I##chi * hd->d##K##F + bd->d##K##chi * hd->d##I##F \
      - bd->gamma##I##K * (                                     \
        bd->gammai11 * bd->d1chi * hd->d1F + bd->gammai22 * bd->d2chi * hd->d2F + bd->gammai33 * bd->d3chi * hd->d3F \
        + bd->gammai12 * bd->d1chi * hd->d2F + bd->gammai13 * bd->d1chi * hd->d3F + bd->gammai23 * bd->d2chi * hd->d3F\
        + bd->gammai21 * bd->d2chi * hd->d1F + bd->gammai31 * bd->d3chi * hd->d1F + bd->gammai32 * bd->d3chi * hd->d2F)) / bd->chi \
  )

#define HORIZON_CALCULATE_DS_TERM3(I, K) \
  (hd->d##I##d##K##F - bd->G##1##I##K * hd->d1F   \
   - bd->G##2##I##K * hd->d2F                     \
   - bd->G##3##I##K * hd->d3F                     \
   - (bd->d##I##chi * hd->d##K##F)  / bd->chi   \
  )


#define HORIZON_CALCULATE_DS1(I, K) \
  (bd->gammai##I##K                 \
   * HORIZON_CALCULATE_DS_TERM3(I, K))


#define HORIZON_CALCULATE_DS0(I, K) \
  (bd->gammai##I##K                 \
   * HORIZON_CALCULATE_DS_TERM2(I, K))

#define HORIZON_CALCULATE_DS(I, K) \
  COSMO_SUMMATION_2_ARGS(HORIZON_CALCULATE_DS_TERM1, I, K)       \
  * HORIZON_CALCULATE_DS_TERM2(I, K)

namespace cosmo
{
 
class Horizon
{
public:

  std::ostream* lstream;

  boost::shared_ptr<tbox::Database>& cosmo_horizon_db;

  const tbox::Dimension& dim;

  Horizon(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const tbox::Dimension& dim_in,
    boost::shared_ptr<tbox::Database> database_in,
    std::ostream* l_stream_in);

  ~Horizon();

  VAR_CREATE(F);

  RK4_IDX_ALL_CREATE(F);
  RK4_PDATA_ALL_CREATE(F);
  RK4_MDA_ACCESS_ALL_CREATE(F);

  VAR_CREATE(s1);
  VAR_CREATE(s2);
  VAR_CREATE(s3);

  RK4_IDX_CREATE(s1, a);
  RK4_PDATA_CREATE(s1, a);
  RK4_MDA_ACCESS_CREATE(s1, a);

  RK4_IDX_CREATE(s2, a);
  RK4_PDATA_CREATE(s2, a);
  RK4_MDA_ACCESS_CREATE(s2, a);

  RK4_IDX_CREATE(s3, a);
  RK4_PDATA_CREATE(s3, a);
  RK4_MDA_ACCESS_CREATE(s3, a);


  void clear(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln);

  
  void addNormVector(
    BSSN *bssn,   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void addNormVector(
    BSSN *bssn, const boost::shared_ptr<hier::PatchLevel> & level);

  void addNormVector(
    BSSN *bssn, const boost::shared_ptr<hier::Patch> & patch);

  
  
  HorizonData getHorizonData(
    idx_t i, idx_t j, idx_t k, BSSNData *bd, const real_t dx[]);  
  
  
  void alloc(const boost::shared_ptr<hier::PatchHierarchy> &hierarchy, idx_t ln);
  
  void addFieldsToList(std::vector<idx_t> &list);

  void registerRKRefinerActive(
    xfer::RefineAlgorithm& refiner,
    boost::shared_ptr<hier::RefineOperator> &space_refine_op);

  
  void registerRKRefiner(
    xfer::RefineAlgorithm& refiner,
    boost::shared_ptr<hier::RefineOperator> &space_refine_op);
  void copyAToP(
    math::HierarchyCellDataOpsReal<real_t> & hcellmath);
  void initPData(
    const boost::shared_ptr<hier::Patch> & patch);
  void initMDA(
    const boost::shared_ptr<hier::Patch> & patch);
  void setLevelTime(
    const boost::shared_ptr<hier::PatchLevel> & level,
    double from_t, double to_t);

  void registerCoarsenActive(
    xfer::CoarsenAlgorithm& coarsener,
    boost::shared_ptr<hier::CoarsenOperator>& coarsen_op);

  
  void K1FinalizePatch(
    const boost::shared_ptr<hier::Patch> & patch);
  void K2FinalizePatch(
    const boost::shared_ptr<hier::Patch> & patch);
  void K3FinalizePatch(
    const boost::shared_ptr<hier::Patch> & patch);
  void K4FinalizePatch(
    const boost::shared_ptr<hier::Patch> & patch);

  void initSurface(const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);

  void stepInit(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);
  
  void RKEvolvePatchBD(const boost::shared_ptr<hier::Patch> & patch, real_t dt);
  real_t getFonBD(
    const boost::shared_ptr<hier::Patch> & patch,
    real_t bx, real_t by, real_t bz,
    hier::Box& g_box);

  bool onTheSurface(idx_t i, idx_t j, idx_t k);

  real_t maxSurfaceMove(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln, idx_t w_idx);


  void RKEvolvePt(
    idx_t i, idx_t j, idx_t k, BSSNData &bd, const real_t dx[], real_t dt);
  void RKEvolveHorizon(
    const boost::shared_ptr<hier::Patch> & patch, BSSN * bssn, real_t dt);

  
  void prepareForK1(
    const boost::shared_ptr<hier::PatchLevel> & level, real_t to_t);
  void prepareForK2(
    const boost::shared_ptr<hier::PatchLevel> & level, real_t to_t);
  void prepareForK3(
    const boost::shared_ptr<hier::PatchLevel> & level, real_t to_t);
  void prepareForK4(
    const boost::shared_ptr<hier::PatchLevel> & level, real_t to_t);

  real_t ev_F(BSSNData *bd, HorizonData *hd, const real_t dx[]);

  void updateBD(
    const boost::shared_ptr<hier::Patch> & patch,
    real_t dt);

  real_t domain_lower[3], domain_upper[3], radius;
  std::vector<double> origin;

  bool is_periodic;
};

} // namespace cosmo


#endif
