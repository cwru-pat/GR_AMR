#ifndef COSMO_GEODESIC
#define COSMO_GEODESIC

#include "../../cosmo_includes.h"
#include "../boundaries/sommerfield.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "../../components/bssn/bssn.h"

#include <SAMRAI/xfer/PatchLevelFillPattern.h>
#include <SAMRAI/xfer/PatchLevelInteriorFillPattern.h>
#include <SAMRAI/xfer/PatchLevelBorderFillPattern.h>
#include <SAMRAI/xfer/PatchLevelEnhancedFillPattern.h>
#include <SAMRAI/xfer/PatchLevelFullFillPattern.h>
#include <SAMRAI/xfer/PatchLevelBorderAndInteriorFillPattern.h>
#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include <SAMRAI/hier/CoarseFineBoundary.h>


#include "geodesic_data.h"
#include "particles.h"

#include "SAMRAI/pdat/IndexVariable.h"
#include "SAMRAI/pdat/IndexData.h"

using namespace SAMRAI;


namespace cosmo
{

#define GEODESIC_DEFINE_CRSPLINES_CHI \
  double a_chi[64], f_chi[64];

#define GEODESIC_DEFINE_CRSPLINES_DCHI(I)       \
  double a_d##I##chi[64], f_d##I##chi[64];
  
#define GEODESIC_DEFINE_CRSPLINES_BETA(I) \
  double a_beta##I[64], f_beta##I[64];

#define GEODESIC_DEFINE_CRSPLINES_DBETA(I,J) \
  double a_d##I##beta##J[64], f_d##I##beta##J[64];
  
#define GEODESIC_DEFINE_CRSPLINES_ALPHA \
  double a_alpha[64], f_alpha[64];

#define GEODESIC_DEFINE_CRSPLINES_DALPHA(I) \
  double a_d##I##alpha[64], f_d##I##alpha[64];

#define GEODESIC_DEFINE_CRSPLINES_M(I,J) \
  double a_m##I##J[64], f_m##I##J[64];

#define GEODESIC_DEFINE_CRSPLINES_DM(K,I,J) \
  double a_d##K##m##I##J[64], f_d##K##m##I##J[64];




#define GEODESIC_CRSPLINES_SET_F_CHI            \
  f_chi[i*16 + j*4 + k] = bd.chi 

#define GEODESIC_CRSPLINES_SET_F_DCHI(I)         \
  f_d##I##chi[i*16 + j*4 + k] = bd.d##I##chi 
  
#define GEODESIC_CRSPLINES_SET_F_BETA(I)        \
  f_beta##I[i*16 + j*4 + k] = bd.beta##I

#define GEODESIC_CRSPLINES_SET_F_DBETA(I,J)             \
  f_d##I##beta##J[i*16 + j*4 + k] = bd.d##I##beta##J

#define GEODESIC_CRSPLINES_SET_F_ALPHA          \
  f_alpha[i*16 + j*4 + k] = bd.DIFFalpha + 1.0 

#define GEODESIC_CRSPLINES_SET_F_DALPHA(I)      \
  f_d##I##alpha[i*16 + j*4 + k] = bd.d##I##a

#define GEODESIC_CRSPLINES_SET_F_M(I,J)         \
  f_m##I##J[i*16 + j*4 + k] = bd.gamma##I##J

#define GEODESIC_CRSPLINES_SET_F_DM(K,I,J)      \
  f_d##K##m##I##J[i*16 + j*4 + k] = bd.d##K##g##I##J


#define GEODESIC_CRSPLINES_CAL_COEF_CHI            \
  compute_tricubic_coeffs(a_chi, f_chi)

#define GEODESIC_CRSPLINES_CAL_COEF_DCHI(I)              \
  compute_tricubic_coeffs(a_d##I##chi, f_d##I##chi)
  
#define GEODESIC_CRSPLINES_CAL_COEF_BETA(I)        \
  compute_tricubic_coeffs(a_beta##I, f_beta##I)

#define GEODESIC_CRSPLINES_CAL_COEF_DBETA(I,J)             \
  compute_tricubic_coeffs(a_d##I##beta##J, f_d##I##beta##J)

#define GEODESIC_CRSPLINES_CAL_COEF_ALPHA          \
  compute_tricubic_coeffs(a_alpha, f_alpha)

#define GEODESIC_CRSPLINES_CAL_COEF_DALPHA(I)      \
  compute_tricubic_coeffs(a_d##I##alpha, f_d##I##alpha)

#define GEODESIC_CRSPLINES_CAL_COEF_M(I,J)         \
  compute_tricubic_coeffs(a_m##I##J, f_m##I##J)

#define GEODESIC_CRSPLINES_CAL_COEF_DM(K,I,J)      \
  compute_tricubic_coeffs(a_d##K##m##I##J, f_d##K##m##I##J)



#define GEODESIC_CRSPLINES_EVAL_CHI            \
  gd->chi = evaluate_interpolation(a_chi, xd, yd, zd)

#define GEODESIC_CRSPLINES_EVAL_DCHI(I)                  \
  gd->d##I##chi = evaluate_interpolation(a_d##I##chi, xd, yd, zd)
  
#define GEODESIC_CRSPLINES_EVAL_BETA(I)        \
  gd->beta##I = evaluate_interpolation(a_beta##I, xd, yd, zd)

#define GEODESIC_CRSPLINES_EVAL_DBETA(I,J)             \
  gd->d##I##beta##J = evaluate_interpolation(a_d##I##beta##J, xd, yd, zd)

#define GEODESIC_CRSPLINES_EVAL_ALPHA          \
  gd->alpha = evaluate_interpolation(a_alpha, xd, yd, zd)

#define GEODESIC_CRSPLINES_EVAL_DALPHA(I)      \
  gd->d##I##alpha = evaluate_interpolation(a_d##I##alpha, xd, yd, zd)

#define GEODESIC_CRSPLINES_EVAL_M(I,J)         \
  gd->m##I##J = evaluate_interpolation(a_m##I##J, xd, yd, zd)

#define GEODESIC_CRSPLINES_EVAL_DM(K,I,J)      \
  gd->d##K##m##I##J = evaluate_interpolation(a_d##K##m##I##J, xd, yd, zd)

#define GEODESIC_CRSPLINES_CAL_DM(K,I,J)      \
  gd->d##K##m##I##J = -2.0 * gd->d##K##chi * gd->m##I##J / pw3(gd->chi) \
    + gd->d##K##m##I##J / pw2(gd->chi)
  


  
class Geodesic
{
 public:
  Geodesic(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const tbox::Dimension& dim_in,
    std::shared_ptr<tbox::Database> database_in,
    std::ostream* l_stream_in,
    real_t KO_damping_coefficient_in,
    int weight_idx);
  
  void initPData(
    const std::shared_ptr<hier::Patch> & patch);

  void insertPatchParticles(
    const std::shared_ptr<hier::Patch> & patch, int src_id, int dst_id);

  void clearParticlesLivingInGhostCells(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int ln);

  void clearParticles(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int idx);
  void clearParticles(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int ln, int idx);
  void clearParticles(
    const std::shared_ptr<hier::Patch> & patch, int idx);

  
  void insertLevelParticles(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int ln, int src_id, int dst_id);

  void preAdvance(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int ln);

  
  void RKEvolvePatch(
    const std::shared_ptr<hier::Patch> & patch, BSSN *bssn, real_t dt);

  void RKEvolveParticle(
    RKParticle &p, GeodesicData &gd, double dt);

  
  void K1FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);

  void K2FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);

  void K3FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);

  void K4FinalizePatch(
    const std::shared_ptr<hier::Patch> & patch);

  void set_gd_values(
  const std::shared_ptr<hier::Patch> & patch, 
  double p_info[], GeodesicData *gd, BSSN *bssn, const real_t dx[], double shift[]);

  
  void clearParticlesCoveredbyFinerLevel(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int ln);

  void particleRedistribution(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    int ln, bool do_update_buffer);

  void printAll(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int idx);

  void geodesic_ic_Schwarzchild_test(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    std::shared_ptr<tbox::Database> cosmo_geodesic_db);

  void geodesic_ic_face_null_test(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    std::shared_ptr<tbox::Database> cosmo_geodesic_db);

  
  void initAll(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy);
  
  void compute_tricubic_coeffs(double *a, double *f);

  double evaluate_interpolation(
    double * a, double x, double y, double z);

  void allocParticles(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int ln);

  void registerRKRefiner(
    xfer::RefineAlgorithm& refiner,
    std::shared_ptr<hier::RefineOperator> &particle_refine_op);

  void regridPreProcessing(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    std::shared_ptr<hier::CoarsenOperator> &particle_coarsen_op);

  void regridPostProcessing(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy);




  std::ostream* lstream;

  const std::shared_ptr<tbox::Database> cosmo_geodesic_db;

  const tbox::Dimension& dim;

  
  std::shared_ptr<pdat::IndexVariable<ParticleContainer,
    pdat::CellGeometry> > pc;

  int pc_idx, pc_s_idx, pc_d_buffer_idx, weight_idx;
  int ghost_width;

  std::shared_ptr<pdat::IndexData<ParticleContainer,
                                  pdat::CellGeometry> > pc_pdata;
  std::shared_ptr<pdat::IndexData<ParticleContainer,
                                  pdat::CellGeometry> > pc_s_pdata;
  std::shared_ptr<pdat::IndexData<ParticleContainer,
                                  pdat::CellGeometry> > pc_d_buffer_pdata;

  double domain_lower[3], domain_upper[3], L[3];
  
  double p0;
  bool save_metric;
  
};
}
#endif