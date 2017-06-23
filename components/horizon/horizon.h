#ifndef COSMO_HORIZON
#define COSMO_HORIZON

#include "../../cosmo_includes.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "horizon_data.h"
#include "../components/bssn/bssn.h"
#include "SAMRAI/hier/LocalId.h"

using namespace SAMRAI;


#define d2d1F d1d2F
#define d3d1F d1d3F
#define d3d2F d2d3F

// christoffel symbols (not conformal)
#define Gc121 Gc112
#define Gc131 Gc113
#define Gc132 Gc123
#define Gc221 Gc212
#define Gc231 Gc213
#define Gc232 Gc223
#define Gc321 Gc312
#define Gc331 Gc313
#define Gc332 Gc323

// christoffel symbols (not conformal) on surface
#define Gs121 Gs112
#define Gs221 Gs212

#define q21 q12

#define mi21 mi12
#define mi31 mi13
#define mi32 mi23

#define R21 R12

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

#define HORIZON_CALCULATE_G(I, J, K) \
  (bd.G##I##J##K - 1.0 / bd.chi * (          \
    (double)(I == K) * bd.d##J##chi                  \
    + (double)(I == J) * bd.d##K##chi                      \
    - (bd.gammai##I##1 * bd.d1chi             \
       + bd.gammai##I##2 * bd.d2chi           \
       + bd.gammai##I##3 * bd.d3chi) * bd.gamma##J##K))

#define HORIZON_DEFINE_TEMP_GI(I, J, K)          \
  double tempGi##I##J##K = 0;

#define HORIZON_DEFINE_TEMP_GJ(I, J, K)          \
  double tempGj##I##J##K = 0;


#define HORIZON_DEFINE_TEMP_MI(I, J)          \
  double tempmi##I##J = 0;

#define HORIZON_DEFINE_TEMP_MJ(I, J)          \
  double tempmj##I##J = 0;


#define HORIZON_DEFINE_TEMP_KI(I, J)          \
  double tempKi##I##J = 0;

#define HORIZON_DEFINE_TEMP_KJ(I, J)          \
  double tempKj##I##J = 0;


#define HORIZON_DEFINE_TEMP_DFI(I)          \
  double tempDFi##I = 0;

#define HORIZON_DEFINE_TEMP_DFJ(I)          \
  double tempDFj##I = 0;


#define HORIZON_DEFINE_TEMP_RI                  \
  double tempRi = 0;

#define HORIZON_DEFINE_TEMP_RJ                  \
  double tempRj = 0;

#define HORIZON_CALCULATE_D1G(I, J, K)          \
  ((G##I##J##K[(theta_i+1)%(2*n_theta)][phi_i] - G##I##J##K[(theta_i + (2*n_theta))% (2*n_theta)][phi_i]) / (PI / n_theta / 2.0)) 

#define HORIZON_CALCULATE_D2G(I, J, K)          \
  ((G##I##J##K[theta_i][(phi_i+1)%(2*n_phi)] - G##I##J##K[theta_i][(phi_i + (2*n_phi))%(2*n_phi)]) / (PI / n_phi)) 


#define HORIZON_INTERPOLATE_G_1(I, J, K)  \
  tempGi##I##J##K += HORIZON_CALCULATE_G(I, J, K) * \
    ( (dx[0] - (x0 - x ) * (2.0 * (i - i0) -1.0)) / dx[0])

#define HORIZON_INTERPOLATE_G_2(I, J, K)  \
  tempGj##I##J##K += tempGi##I##J##K * \
    ( (dx[1] - (y0 - y ) * (2.0 * (j - j0) -1.0))/ dx[1])

#define HORIZON_INTERPOLATE_G_3(I, J, K)  \
  kd->Gc##I##J##K += tempGj##I##J##K * \
    ((dx[2] - (z0 - z ) * (2.0 * (k - k0) -1.0)) / dx[2])

#define HORIZON_INTERPOLATE_M_1(I, J)  \
  tempmi##I##J += 1.0 / PW2(bd.chi) * bd.gamma##I##J *        \
    ( (dx[0] - (x0 - x ) * (2.0 * (i - i0) -1.0)) / dx[0])

#define HORIZON_INTERPOLATE_M_2(I, J)  \
  tempmj##I##J += tempmi##I##J * \
    ( (dx[1] - (y0 - y ) * (2.0 * (j - j0) -1.0))/ dx[1])

#define HORIZON_INTERPOLATE_M_3(I, J)  \
  kd->m##I##J += tempmj##I##J * \
    ((dx[2] - (z0 - z ) * (2.0 * (k - k0) -1.0)) / dx[2])

#define HORIZON_INTERPOLATE_K_1(I, J)  \
  tempKi##I##J += 1.0 / PW2(bd.chi) * (bd.A##I##J +  bd.gamma##I##J * bd.K / 3.0) * \
    ( (dx[0] - (x0 - x ) * (2.0 * (i - i0) -1.0)) / dx[0])

#define HORIZON_INTERPOLATE_K_2(I, J)  \
  tempKj##I##J += tempKi##I##J * \
    ( (dx[1] - (y0 - y ) * (2.0 * (j - j0) -1.0))/ dx[1])


#define HORIZON_INTERPOLATE_K_3(I, J)  \
  kd->K##I##J += tempKj##I##J * \
    ((dx[2] - (z0 - z ) * (2.0 * (k - k0) -1.0)) / dx[2])



#define HORIZON_INTERPOLATE_DF_1(I)  \
  tempDFi##I += hd.d##I##F * \
    ( (dx[0] - (x0 - x ) * (2.0 * (i - i0) -1.0)) / dx[0])

#define HORIZON_INTERPOLATE_DF_2(I)  \
  tempDFj##I += tempDFi##I * \
    ( (dx[1] - (y0 - y ) * (2.0 * (j - j0) -1.0))/ dx[1])


#define HORIZON_INTERPOLATE_DF_3(I)  \
  kd->d##I##F += tempDFj##I * \
    ((dx[2] - (z0 - z ) * (2.0 * (k - k0) -1.0)) / dx[2])



#define HORIZON_INTERPOLATE_R_1                 \
  tempRi += bd.ricci *                  \
    ( (dx[0] - (x0 - x ) * (2.0 * (i - i0) -1.0)) / dx[0])    

#define HORIZON_INTERPOLATE_R_2                 \
  tempRj += tempRi *                          \
    ( (dx[1] - (y0 - y ) * (2.0 * (j - j0) -1.0))/ dx[1])

#define HORIZON_INTERPOLATE_R_3                         \
  kd->R += tempRj *                                     \
    ((dx[2] - (z0 - z ) * (2.0 * (k - k0) -1.0)) / dx[2])




namespace cosmo
{


struct Cell{
  int i, j, k;
  int patch_id;
  int level_num;  
  double radius;
  double d; // distance on the surface from the certain direction
};
 
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
    std::ostream* l_stream_in,
    int w_idx_in);

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
  bool initSphericalSurface(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn,   boost::shared_ptr<hier::RefineOperator> space_refine_op);

  void stepInit(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);
  
  void RKEvolvePatchBD(const boost::shared_ptr<hier::Patch> & patch, real_t dt);
  real_t getFonBD(
    const boost::shared_ptr<hier::Patch> & patch,
    real_t bx, real_t by, real_t bz,
    hier::Box& g_box);

  bool onTheSurface(idx_t i, idx_t j, idx_t k);
  bool belowTheSurface(idx_t i, idx_t j, idx_t k);
  
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

  void findKilling(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn);
  void initGridding(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy);
  
  void initG(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn);

  real_t findMaxHorizonRadius(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    double theta_0, double phi_0);
  void findM(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,  double x[], BSSN *bssn);
  void transportKillingPhi(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    idx_t theta_i, idx_t phi_f, double k_theta_0, double k_phi_0, double k_L_0, BSSN * bssn);
  void transportKillingTheta(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    idx_t phi_i, double k_theta_0, double k_phi_0, double k_L_0, BSSN * bssn);

  real_t findRadius(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    double theta_0, double phi_0);
  void findPatch(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    double theta_0, double phi_0);

  real_t getRadius(int theta_i, int phi_i);
  
  void set_kd_values(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    double theta, double phi, int theta_i, int phi_i, double r, KillingData *kd, BSSN * bssn);

  void set_G_values(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    double theta, double phi, int theta_i, int phi_i, double r, KillingData *kd, BSSN * bssn);

  void set_norm_values(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    double theta, double phi, int theta_i, int phi_i, double r, KillingData *kd, BSSN * bssn);
  real_t angularMomentum(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn);
  real_t area(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn);

  
  real_t ev_k_theta_dtheta(KillingData *kd, int theta_i, int phi_i);
  real_t ev_k_phi_dtheta(KillingData *kd, int theta_i, int phi_i);
  real_t ev_k_L_dtheta(KillingData *kd, int theta_i, int phi_i);
  real_t ev_k_theta_dphi(KillingData *kd, int theta_i, int phi_i);
  real_t ev_k_phi_dphi(KillingData *kd, int theta_i, int phi_i);
  real_t ev_k_L_dphi(KillingData *kd, int theta_i, int phi_i);

  void convertToVector(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn);
  real_t getNormFactor();
  void normKilling();
  real_t interp_k_phi(
    double theta, double phi);
  real_t interp_k_theta(
    double theta, double phi);





  
  real_t domain_lower[3], domain_upper[3], radius;
  std::vector<double> origin;

  // whether or not periodic boundary is used,
  // if so, will initialize 2 initial surface
  bool is_periodic;

  bool is_sphere;

  double const_radius;
  
  double radius_limit;
  
  std::vector<std::vector<double>> k_theta, k_phi, k_L;
  idx_t n_theta, n_phi;

  int patch_work_i, patch_work_j, patch_work_k, patch_work_level, local_id;
  hier::GlobalId patch_work_id;
  idx_t patch_work_mpi_rank, cur_mpi_rank;
  real_t min_d;
  double M[3][3];
  idx_t w_idx;
  hier::GlobalId invalid_id;
  std::vector<std::vector<double>> G111, G112, G122, G211, G212, G222, ah_radius;

};

} // namespace cosmo


#endif
