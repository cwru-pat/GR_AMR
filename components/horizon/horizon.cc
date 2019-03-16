#include "horizon.h"
#include "../../cosmo_includes.h"
#include "../../utils/math.h"
#include "../bssn/bssn_data.h"
#include <complex.h>
#include "../../utils/Eigen/Dense"

using namespace SAMRAI;

namespace cosmo
{
  
HorizonStatistics::HorizonStatistics(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  const tbox::Dimension& dim_in,
  std::shared_ptr<tbox::Database> database_in,
  int w_idx_in,
  AHFinderDirect::Horizon *horizon_in):
  cosmo_horizon_db(database_in),
  dim(dim_in),
  n_phi(cosmo_horizon_db->getIntegerWithDefault("n_phi", 36)),
  w_idx(w_idx_in),
  invalid_id( hier::LocalId::getInvalidId(), tbox::SAMRAI_MPI::getInvalidRank()),
  non_zero_angular_momentum(cosmo_horizon_db->getBoolWithDefault("non_zero_angular_momentum", true)),
  horizon(horizon_in)
{
  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;


  const double * lower = &grid_geometry.getXLower()[0];
  const double * upper = &grid_geometry.getXUpper()[0];
  for(int i = 0 ; i < DIM; i++)
  {
    domain_lower[i] = lower[i];
    domain_upper[i] = upper[i];
  }

  origin.resize(3); // position of the locale of origin of horizon
  coord_origin.resize(3); // position of the coordinate origin where spherical coordinate is built

  for(int i = 0; i < 3; i ++)
    coord_origin[i] = (upper[i] - lower[i]) / 2.0;
  
  patch_work_i = patch_work_j = patch_work_k = -1;
  
}
HorizonStatistics::~HorizonStatistics()
{
  
}

void HorizonStatistics::compute_tricubic_coeffs(double *a, double *f)
{
  a[0] = f[21];
  a[1] = -f[20]/3. - f[21]/2. + f[22] - f[23]/6.;
  a[2] = f[20]/2. - f[21] + f[22]/2.;
  a[3] = -f[20]/6. + f[21]/2. - f[22]/2. + f[23]/6.;
  a[4] = -f[17]/3. - f[21]/2. + f[25] - f[29]/6.;
  a[5] = f[16]/9. + f[17]/6. - f[18]/3. + f[19]/18. + f[20]/6. + f[21]/4. - f[22]/2. + f[23]/12. - f[24]/3. - f[25]/2. + f[26] - f[27]/6. + f[28]/18. + f[29]/12. - f[30]/6. + f[31]/36.;
  a[6] = -f[16]/6. + f[17]/3. - f[18]/6. - f[20]/4. + f[21]/2. - f[22]/4. + f[24]/2. - f[25] + f[26]/2. - f[28]/12. + f[29]/6. - f[30]/12.;
  a[7] = f[16]/18. - f[17]/6. + f[18]/6. - f[19]/18. + f[20]/12. - f[21]/4. + f[22]/4. - f[23]/12. - f[24]/6. + f[25]/2. - f[26]/2. + f[27]/6. + f[28]/36. - f[29]/12. + f[30]/12. - f[31]/36.;
  a[8] = f[17]/2. - f[21] + f[25]/2.;
  a[9] = -f[16]/6. - f[17]/4. + f[18]/2. - f[19]/12. + f[20]/3. + f[21]/2. - f[22] + f[23]/6. - f[24]/6. - f[25]/4. + f[26]/2. - f[27]/12.;
  a[10] = f[16]/4. - f[17]/2. + f[18]/4. - f[20]/2. + f[21] - f[22]/2. + f[24]/4. - f[25]/2. + f[26]/4.;
  a[11] = -f[16]/12. + f[17]/4. - f[18]/4. + f[19]/12. + f[20]/6. - f[21]/2. + f[22]/2. - f[23]/6. - f[24]/12. + f[25]/4. - f[26]/4. + f[27]/12.;
  a[12] = -f[17]/6. + f[21]/2. - f[25]/2. + f[29]/6.;
  a[13] = f[16]/18. + f[17]/12. - f[18]/6. + f[19]/36. - f[20]/6. - f[21]/4. + f[22]/2. - f[23]/12. + f[24]/6. + f[25]/4. - f[26]/2. + f[27]/12. - f[28]/18. - f[29]/12. + f[30]/6. - f[31]/36.;
  a[14] = -f[16]/12. + f[17]/6. - f[18]/12. + f[20]/4. - f[21]/2. + f[22]/4. - f[24]/4. + f[25]/2. - f[26]/4. + f[28]/12. - f[29]/6. + f[30]/12.;
  a[15] = f[16]/36. - f[17]/12. + f[18]/12. - f[19]/36. - f[20]/12. + f[21]/4. - f[22]/4. + f[23]/12. + f[24]/12. - f[25]/4. + f[26]/4. - f[27]/12. - f[28]/36. + f[29]/12. - f[30]/12. + f[31]/36.;
  a[16] = -f[5]/3. - f[21]/2. + f[37] - f[53]/6.;
  a[17] = f[4]/9. + f[5]/6. - f[6]/3. + f[7]/18. + f[20]/6. + f[21]/4. - f[22]/2. + f[23]/12. - f[36]/3. - f[37]/2. + f[38] - f[39]/6. + f[52]/18. + f[53]/12. - f[54]/6. + f[55]/36.;
  a[18] = -f[4]/6. + f[5]/3. - f[6]/6. - f[20]/4. + f[21]/2. - f[22]/4. + f[36]/2. - f[37] + f[38]/2. - f[52]/12. + f[53]/6. - f[54]/12.;
  a[19] = f[4]/18. - f[5]/6. + f[6]/6. - f[7]/18. + f[20]/12. - f[21]/4. + f[22]/4. - f[23]/12. - f[36]/6. + f[37]/2. - f[38]/2. + f[39]/6. + f[52]/36. - f[53]/12. + f[54]/12. - f[55]/36.;
  a[20] = f[1]/9. + f[5]/6. - f[9]/3. + f[13]/18. + f[17]/6. + f[21]/4. - f[25]/2. + f[29]/12. - f[33]/3. - f[37]/2. + f[41] - f[45]/6. + f[49]/18. + f[53]/12. - f[57]/6. + f[61]/36.;
  a[21] = -f[0]/27. - f[1]/18. + f[2]/9. - f[3]/54. - f[4]/18. - f[5]/12. + f[6]/6. - f[7]/36. + f[8]/9. + f[9]/6. - f[10]/3. + f[11]/18. - f[12]/54. - f[13]/36. + f[14]/18. - f[15]/108. - f[16]/18. - f[17]/12. + f[18]/6. - f[19]/36. - f[20]/12. - f[21]/8. + f[22]/4. - f[23]/24. + f[24]/6. + f[25]/4. - f[26]/2. + f[27]/12. - f[28]/36. - f[29]/24. + f[30]/12. - f[31]/72. + f[32]/9. + f[33]/6. - f[34]/3. + f[35]/18. + f[36]/6. + f[37]/4. - f[38]/2. + f[39]/12. - f[40]/3. - f[41]/2. + f[42] - f[43]/6. + f[44]/18. + f[45]/12. - f[46]/6. + f[47]/36. - f[48]/54. - f[49]/36. + f[50]/18. - f[51]/108. - f[52]/36. - f[53]/24. + f[54]/12. - f[55]/72. + f[56]/18. + f[57]/12. - f[58]/6. + f[59]/36. - f[60]/108. - f[61]/72. + f[62]/36. - f[63]/216.;
  a[22] = f[0]/18. - f[1]/9. + f[2]/18. + f[4]/12. - f[5]/6. + f[6]/12. - f[8]/6. + f[9]/3. - f[10]/6. + f[12]/36. - f[13]/18. + f[14]/36. + f[16]/12. - f[17]/6. + f[18]/12. + f[20]/8. - f[21]/4. + f[22]/8. - f[24]/4. + f[25]/2. - f[26]/4. + f[28]/24. - f[29]/12. + f[30]/24. - f[32]/6. + f[33]/3. - f[34]/6. - f[36]/4. + f[37]/2. - f[38]/4. + f[40]/2. - f[41] + f[42]/2. - f[44]/12. + f[45]/6. - f[46]/12. + f[48]/36. - f[49]/18. + f[50]/36. + f[52]/24. - f[53]/12. + f[54]/24. - f[56]/12. + f[57]/6. - f[58]/12. + f[60]/72. - f[61]/36. + f[62]/72.;
  a[23] = -f[0]/54. + f[1]/18. - f[2]/18. + f[3]/54. - f[4]/36. + f[5]/12. - f[6]/12. + f[7]/36. + f[8]/18. - f[9]/6. + f[10]/6. - f[11]/18. - f[12]/108. + f[13]/36. - f[14]/36. + f[15]/108. - f[16]/36. + f[17]/12. - f[18]/12. + f[19]/36. - f[20]/24. + f[21]/8. - f[22]/8. + f[23]/24. + f[24]/12. - f[25]/4. + f[26]/4. - f[27]/12. - f[28]/72. + f[29]/24. - f[30]/24. + f[31]/72. + f[32]/18. - f[33]/6. + f[34]/6. - f[35]/18. + f[36]/12. - f[37]/4. + f[38]/4. - f[39]/12. - f[40]/6. + f[41]/2. - f[42]/2. + f[43]/6. + f[44]/36. - f[45]/12. + f[46]/12. - f[47]/36. - f[48]/108. + f[49]/36. - f[50]/36. + f[51]/108. - f[52]/72. + f[53]/24. - f[54]/24. + f[55]/72. + f[56]/36. - f[57]/12. + f[58]/12. - f[59]/36. - f[60]/216. + f[61]/72. - f[62]/72. + f[63]/216.;
  a[24] = -f[1]/6. + f[5]/3. - f[9]/6. - f[17]/4. + f[21]/2. - f[25]/4. + f[33]/2. - f[37] + f[41]/2. - f[49]/12. + f[53]/6. - f[57]/12.;
  a[25] = f[0]/18. + f[1]/12. - f[2]/6. + f[3]/36. - f[4]/9. - f[5]/6. + f[6]/3. - f[7]/18. + f[8]/18. + f[9]/12. - f[10]/6. + f[11]/36. + f[16]/12. + f[17]/8. - f[18]/4. + f[19]/24. - f[20]/6. - f[21]/4. + f[22]/2. - f[23]/12. + f[24]/12. + f[25]/8. - f[26]/4. + f[27]/24. - f[32]/6. - f[33]/4. + f[34]/2. - f[35]/12. + f[36]/3. + f[37]/2. - f[38] + f[39]/6. - f[40]/6. - f[41]/4. + f[42]/2. - f[43]/12. + f[48]/36. + f[49]/24. - f[50]/12. + f[51]/72. - f[52]/18. - f[53]/12. + f[54]/6. - f[55]/36. + f[56]/36. + f[57]/24. - f[58]/12. + f[59]/72.;
  a[26] = -f[0]/12. + f[1]/6. - f[2]/12. + f[4]/6. - f[5]/3. + f[6]/6. - f[8]/12. + f[9]/6. - f[10]/12. - f[16]/8. + f[17]/4. - f[18]/8. + f[20]/4. - f[21]/2. + f[22]/4. - f[24]/8. + f[25]/4. - f[26]/8. + f[32]/4. - f[33]/2. + f[34]/4. - f[36]/2. + f[37] - f[38]/2. + f[40]/4. - f[41]/2. + f[42]/4. - f[48]/24. + f[49]/12. - f[50]/24. + f[52]/12. - f[53]/6. + f[54]/12. - f[56]/24. + f[57]/12. - f[58]/24.;
  a[27] = f[0]/36. - f[1]/12. + f[2]/12. - f[3]/36. - f[4]/18. + f[5]/6. - f[6]/6. + f[7]/18. + f[8]/36. - f[9]/12. + f[10]/12. - f[11]/36. + f[16]/24. - f[17]/8. + f[18]/8. - f[19]/24. - f[20]/12. + f[21]/4. - f[22]/4. + f[23]/12. + f[24]/24. - f[25]/8. + f[26]/8. - f[27]/24. - f[32]/12. + f[33]/4. - f[34]/4. + f[35]/12. + f[36]/6. - f[37]/2. + f[38]/2. - f[39]/6. - f[40]/12. + f[41]/4. - f[42]/4. + f[43]/12. + f[48]/72. - f[49]/24. + f[50]/24. - f[51]/72. - f[52]/36. + f[53]/12. - f[54]/12. + f[55]/36. + f[56]/72. - f[57]/24. + f[58]/24. - f[59]/72.;
  a[28] = f[1]/18. - f[5]/6. + f[9]/6. - f[13]/18. + f[17]/12. - f[21]/4. + f[25]/4. - f[29]/12. - f[33]/6. + f[37]/2. - f[41]/2. + f[45]/6. + f[49]/36. - f[53]/12. + f[57]/12. - f[61]/36.;
  a[29] = -f[0]/54. - f[1]/36. + f[2]/18. - f[3]/108. + f[4]/18. + f[5]/12. - f[6]/6. + f[7]/36. - f[8]/18. - f[9]/12. + f[10]/6. - f[11]/36. + f[12]/54. + f[13]/36. - f[14]/18. + f[15]/108. - f[16]/36. - f[17]/24. + f[18]/12. - f[19]/72. + f[20]/12. + f[21]/8. - f[22]/4. + f[23]/24. - f[24]/12. - f[25]/8. + f[26]/4. - f[27]/24. + f[28]/36. + f[29]/24. - f[30]/12. + f[31]/72. + f[32]/18. + f[33]/12. - f[34]/6. + f[35]/36. - f[36]/6. - f[37]/4. + f[38]/2. - f[39]/12. + f[40]/6. + f[41]/4. - f[42]/2. + f[43]/12. - f[44]/18. - f[45]/12. + f[46]/6. - f[47]/36. - f[48]/108. - f[49]/72. + f[50]/36. - f[51]/216. + f[52]/36. + f[53]/24. - f[54]/12. + f[55]/72. - f[56]/36. - f[57]/24. + f[58]/12. - f[59]/72. + f[60]/108. + f[61]/72. - f[62]/36. + f[63]/216.;
  a[30] = f[0]/36. - f[1]/18. + f[2]/36. - f[4]/12. + f[5]/6. - f[6]/12. + f[8]/12. - f[9]/6. + f[10]/12. - f[12]/36. + f[13]/18. - f[14]/36. + f[16]/24. - f[17]/12. + f[18]/24. - f[20]/8. + f[21]/4. - f[22]/8. + f[24]/8. - f[25]/4. + f[26]/8. - f[28]/24. + f[29]/12. - f[30]/24. - f[32]/12. + f[33]/6. - f[34]/12. + f[36]/4. - f[37]/2. + f[38]/4. - f[40]/4. + f[41]/2. - f[42]/4. + f[44]/12. - f[45]/6. + f[46]/12. + f[48]/72. - f[49]/36. + f[50]/72. - f[52]/24. + f[53]/12. - f[54]/24. + f[56]/24. - f[57]/12. + f[58]/24. - f[60]/72. + f[61]/36. - f[62]/72.;
  a[31] = -f[0]/108. + f[1]/36. - f[2]/36. + f[3]/108. + f[4]/36. - f[5]/12. + f[6]/12. - f[7]/36. - f[8]/36. + f[9]/12. - f[10]/12. + f[11]/36. + f[12]/108. - f[13]/36. + f[14]/36. - f[15]/108. - f[16]/72. + f[17]/24. - f[18]/24. + f[19]/72. + f[20]/24. - f[21]/8. + f[22]/8. - f[23]/24. - f[24]/24. + f[25]/8. - f[26]/8. + f[27]/24. + f[28]/72. - f[29]/24. + f[30]/24. - f[31]/72. + f[32]/36. - f[33]/12. + f[34]/12. - f[35]/36. - f[36]/12. + f[37]/4. - f[38]/4. + f[39]/12. + f[40]/12. - f[41]/4. + f[42]/4. - f[43]/12. - f[44]/36. + f[45]/12. - f[46]/12. + f[47]/36. - f[48]/216. + f[49]/72. - f[50]/72. + f[51]/216. + f[52]/72. - f[53]/24. + f[54]/24. - f[55]/72. - f[56]/72. + f[57]/24. - f[58]/24. + f[59]/72. + f[60]/216. - f[61]/72. + f[62]/72. - f[63]/216.;
  a[32] = f[5]/2. - f[21] + f[37]/2.;
  a[33] = -f[4]/6. - f[5]/4. + f[6]/2. - f[7]/12. + f[20]/3. + f[21]/2. - f[22] + f[23]/6. - f[36]/6. - f[37]/4. + f[38]/2. - f[39]/12.;
  a[34] = f[4]/4. - f[5]/2. + f[6]/4. - f[20]/2. + f[21] - f[22]/2. + f[36]/4. - f[37]/2. + f[38]/4.;
  a[35] = -f[4]/12. + f[5]/4. - f[6]/4. + f[7]/12. + f[20]/6. - f[21]/2. + f[22]/2. - f[23]/6. - f[36]/12. + f[37]/4. - f[38]/4. + f[39]/12.;
  a[36] = -f[1]/6. - f[5]/4. + f[9]/2. - f[13]/12. + f[17]/3. + f[21]/2. - f[25] + f[29]/6. - f[33]/6. - f[37]/4. + f[41]/2. - f[45]/12.;
  a[37] = f[0]/18. + f[1]/12. - f[2]/6. + f[3]/36. + f[4]/12. + f[5]/8. - f[6]/4. + f[7]/24. - f[8]/6. - f[9]/4. + f[10]/2. - f[11]/12. + f[12]/36. + f[13]/24. - f[14]/12. + f[15]/72. - f[16]/9. - f[17]/6. + f[18]/3. - f[19]/18. - f[20]/6. - f[21]/4. + f[22]/2. - f[23]/12. + f[24]/3. + f[25]/2. - f[26] + f[27]/6. - f[28]/18. - f[29]/12. + f[30]/6. - f[31]/36. + f[32]/18. + f[33]/12. - f[34]/6. + f[35]/36. + f[36]/12. + f[37]/8. - f[38]/4. + f[39]/24. - f[40]/6. - f[41]/4. + f[42]/2. - f[43]/12. + f[44]/36. + f[45]/24. - f[46]/12. + f[47]/72.;
  a[38] = -f[0]/12. + f[1]/6. - f[2]/12. - f[4]/8. + f[5]/4. - f[6]/8. + f[8]/4. - f[9]/2. + f[10]/4. - f[12]/24. + f[13]/12. - f[14]/24. + f[16]/6. - f[17]/3. + f[18]/6. + f[20]/4. - f[21]/2. + f[22]/4. - f[24]/2. + f[25] - f[26]/2. + f[28]/12. - f[29]/6. + f[30]/12. - f[32]/12. + f[33]/6. - f[34]/12. - f[36]/8. + f[37]/4. - f[38]/8. + f[40]/4. - f[41]/2. + f[42]/4. - f[44]/24. + f[45]/12. - f[46]/24.;
  a[39] = f[0]/36. - f[1]/12. + f[2]/12. - f[3]/36. + f[4]/24. - f[5]/8. + f[6]/8. - f[7]/24. - f[8]/12. + f[9]/4. - f[10]/4. + f[11]/12. + f[12]/72. - f[13]/24. + f[14]/24. - f[15]/72. - f[16]/18. + f[17]/6. - f[18]/6. + f[19]/18. - f[20]/12. + f[21]/4. - f[22]/4. + f[23]/12. + f[24]/6. - f[25]/2. + f[26]/2. - f[27]/6. - f[28]/36. + f[29]/12. - f[30]/12. + f[31]/36. + f[32]/36. - f[33]/12. + f[34]/12. - f[35]/36. + f[36]/24. - f[37]/8. + f[38]/8. - f[39]/24. - f[40]/12. + f[41]/4. - f[42]/4. + f[43]/12. + f[44]/72. - f[45]/24. + f[46]/24. - f[47]/72.;
  a[40] = f[1]/4. - f[5]/2. + f[9]/4. - f[17]/2. + f[21] - f[25]/2. + f[33]/4. - f[37]/2. + f[41]/4.;
  a[41] = -f[0]/12. - f[1]/8. + f[2]/4. - f[3]/24. + f[4]/6. + f[5]/4. - f[6]/2. + f[7]/12. - f[8]/12. - f[9]/8. + f[10]/4. - f[11]/24. + f[16]/6. + f[17]/4. - f[18]/2. + f[19]/12. - f[20]/3. - f[21]/2. + f[22] - f[23]/6. + f[24]/6. + f[25]/4. - f[26]/2. + f[27]/12. - f[32]/12. - f[33]/8. + f[34]/4. - f[35]/24. + f[36]/6. + f[37]/4. - f[38]/2. + f[39]/12. - f[40]/12. - f[41]/8. + f[42]/4. - f[43]/24.;
  a[42] = f[0]/8. - f[1]/4. + f[2]/8. - f[4]/4. + f[5]/2. - f[6]/4. + f[8]/8. - f[9]/4. + f[10]/8. - f[16]/4. + f[17]/2. - f[18]/4. + f[20]/2. - f[21] + f[22]/2. - f[24]/4. + f[25]/2. - f[26]/4. + f[32]/8. - f[33]/4. + f[34]/8. - f[36]/4. + f[37]/2. - f[38]/4. + f[40]/8. - f[41]/4. + f[42]/8.;
  a[43] = -f[0]/24. + f[1]/8. - f[2]/8. + f[3]/24. + f[4]/12. - f[5]/4. + f[6]/4. - f[7]/12. - f[8]/24. + f[9]/8. - f[10]/8. + f[11]/24. + f[16]/12. - f[17]/4. + f[18]/4. - f[19]/12. - f[20]/6. + f[21]/2. - f[22]/2. + f[23]/6. + f[24]/12. - f[25]/4. + f[26]/4. - f[27]/12. - f[32]/24. + f[33]/8. - f[34]/8. + f[35]/24. + f[36]/12. - f[37]/4. + f[38]/4. - f[39]/12. - f[40]/24. + f[41]/8. - f[42]/8. + f[43]/24.;
  a[44] = -f[1]/12. + f[5]/4. - f[9]/4. + f[13]/12. + f[17]/6. - f[21]/2. + f[25]/2. - f[29]/6. - f[33]/12. + f[37]/4. - f[41]/4. + f[45]/12.;
  a[45] = f[0]/36. + f[1]/24. - f[2]/12. + f[3]/72. - f[4]/12. - f[5]/8. + f[6]/4. - f[7]/24. + f[8]/12. + f[9]/8. - f[10]/4. + f[11]/24. - f[12]/36. - f[13]/24. + f[14]/12. - f[15]/72. - f[16]/18. - f[17]/12. + f[18]/6. - f[19]/36. + f[20]/6. + f[21]/4. - f[22]/2. + f[23]/12. - f[24]/6. - f[25]/4. + f[26]/2. - f[27]/12. + f[28]/18. + f[29]/12. - f[30]/6. + f[31]/36. + f[32]/36. + f[33]/24. - f[34]/12. + f[35]/72. - f[36]/12. - f[37]/8. + f[38]/4. - f[39]/24. + f[40]/12. + f[41]/8. - f[42]/4. + f[43]/24. - f[44]/36. - f[45]/24. + f[46]/12. - f[47]/72.;
  a[46] = -f[0]/24. + f[1]/12. - f[2]/24. + f[4]/8. - f[5]/4. + f[6]/8. - f[8]/8. + f[9]/4. - f[10]/8. + f[12]/24. - f[13]/12. + f[14]/24. + f[16]/12. - f[17]/6. + f[18]/12. - f[20]/4. + f[21]/2. - f[22]/4. + f[24]/4. - f[25]/2. + f[26]/4. - f[28]/12. + f[29]/6. - f[30]/12. - f[32]/24. + f[33]/12. - f[34]/24. + f[36]/8. - f[37]/4. + f[38]/8. - f[40]/8. + f[41]/4. - f[42]/8. + f[44]/24. - f[45]/12. + f[46]/24.;
  a[47] = f[0]/72. - f[1]/24. + f[2]/24. - f[3]/72. - f[4]/24. + f[5]/8. - f[6]/8. + f[7]/24. + f[8]/24. - f[9]/8. + f[10]/8. - f[11]/24. - f[12]/72. + f[13]/24. - f[14]/24. + f[15]/72. - f[16]/36. + f[17]/12. - f[18]/12. + f[19]/36. + f[20]/12. - f[21]/4. + f[22]/4. - f[23]/12. - f[24]/12. + f[25]/4. - f[26]/4. + f[27]/12. + f[28]/36. - f[29]/12. + f[30]/12. - f[31]/36. + f[32]/72. - f[33]/24. + f[34]/24. - f[35]/72. - f[36]/24. + f[37]/8. - f[38]/8. + f[39]/24. + f[40]/24. - f[41]/8. + f[42]/8. - f[43]/24. - f[44]/72. + f[45]/24. - f[46]/24. + f[47]/72.;
  a[48] = -f[5]/6. + f[21]/2. - f[37]/2. + f[53]/6.;
  a[49] = f[4]/18. + f[5]/12. - f[6]/6. + f[7]/36. - f[20]/6. - f[21]/4. + f[22]/2. - f[23]/12. + f[36]/6. + f[37]/4. - f[38]/2. + f[39]/12. - f[52]/18. - f[53]/12. + f[54]/6. - f[55]/36.;
  a[50] = -f[4]/12. + f[5]/6. - f[6]/12. + f[20]/4. - f[21]/2. + f[22]/4. - f[36]/4. + f[37]/2. - f[38]/4. + f[52]/12. - f[53]/6. + f[54]/12.;
  a[51] = f[4]/36. - f[5]/12. + f[6]/12. - f[7]/36. - f[20]/12. + f[21]/4. - f[22]/4. + f[23]/12. + f[36]/12. - f[37]/4. + f[38]/4. - f[39]/12. - f[52]/36. + f[53]/12. - f[54]/12. + f[55]/36.;
  a[52] = f[1]/18. + f[5]/12. - f[9]/6. + f[13]/36. - f[17]/6. - f[21]/4. + f[25]/2. - f[29]/12. + f[33]/6. + f[37]/4. - f[41]/2. + f[45]/12. - f[49]/18. - f[53]/12. + f[57]/6. - f[61]/36.;
  a[53] = -f[0]/54. - f[1]/36. + f[2]/18. - f[3]/108. - f[4]/36. - f[5]/24. + f[6]/12. - f[7]/72. + f[8]/18. + f[9]/12. - f[10]/6. + f[11]/36. - f[12]/108. - f[13]/72. + f[14]/36. - f[15]/216. + f[16]/18. + f[17]/12. - f[18]/6. + f[19]/36. + f[20]/12. + f[21]/8. - f[22]/4. + f[23]/24. - f[24]/6. - f[25]/4. + f[26]/2. - f[27]/12. + f[28]/36. + f[29]/24. - f[30]/12. + f[31]/72. - f[32]/18. - f[33]/12. + f[34]/6. - f[35]/36. - f[36]/12. - f[37]/8. + f[38]/4. - f[39]/24. + f[40]/6. + f[41]/4. - f[42]/2. + f[43]/12. - f[44]/36. - f[45]/24. + f[46]/12. - f[47]/72. + f[48]/54. + f[49]/36. - f[50]/18. + f[51]/108. + f[52]/36. + f[53]/24. - f[54]/12. + f[55]/72. - f[56]/18. - f[57]/12. + f[58]/6. - f[59]/36. + f[60]/108. + f[61]/72. - f[62]/36. + f[63]/216.;
  a[54] = f[0]/36. - f[1]/18. + f[2]/36. + f[4]/24. - f[5]/12. + f[6]/24. - f[8]/12. + f[9]/6. - f[10]/12. + f[12]/72. - f[13]/36. + f[14]/72. - f[16]/12. + f[17]/6. - f[18]/12. - f[20]/8. + f[21]/4. - f[22]/8. + f[24]/4. - f[25]/2. + f[26]/4. - f[28]/24. + f[29]/12. - f[30]/24. + f[32]/12. - f[33]/6. + f[34]/12. + f[36]/8. - f[37]/4. + f[38]/8. - f[40]/4. + f[41]/2. - f[42]/4. + f[44]/24. - f[45]/12. + f[46]/24. - f[48]/36. + f[49]/18. - f[50]/36. - f[52]/24. + f[53]/12. - f[54]/24. + f[56]/12. - f[57]/6. + f[58]/12. - f[60]/72. + f[61]/36. - f[62]/72.;
  a[55] = -f[0]/108. + f[1]/36. - f[2]/36. + f[3]/108. - f[4]/72. + f[5]/24. - f[6]/24. + f[7]/72. + f[8]/36. - f[9]/12. + f[10]/12. - f[11]/36. - f[12]/216. + f[13]/72. - f[14]/72. + f[15]/216. + f[16]/36. - f[17]/12. + f[18]/12. - f[19]/36. + f[20]/24. - f[21]/8. + f[22]/8. - f[23]/24. - f[24]/12. + f[25]/4. - f[26]/4. + f[27]/12. + f[28]/72. - f[29]/24. + f[30]/24. - f[31]/72. - f[32]/36. + f[33]/12. - f[34]/12. + f[35]/36. - f[36]/24. + f[37]/8. - f[38]/8. + f[39]/24. + f[40]/12. - f[41]/4. + f[42]/4. - f[43]/12. - f[44]/72. + f[45]/24. - f[46]/24. + f[47]/72. + f[48]/108. - f[49]/36. + f[50]/36. - f[51]/108. + f[52]/72. - f[53]/24. + f[54]/24. - f[55]/72. - f[56]/36. + f[57]/12. - f[58]/12. + f[59]/36. + f[60]/216. - f[61]/72. + f[62]/72. - f[63]/216.;
  a[56] = -f[1]/12. + f[5]/6. - f[9]/12. + f[17]/4. - f[21]/2. + f[25]/4. - f[33]/4. + f[37]/2. - f[41]/4. + f[49]/12. - f[53]/6. + f[57]/12.;
  a[57] = f[0]/36. + f[1]/24. - f[2]/12. + f[3]/72. - f[4]/18. - f[5]/12. + f[6]/6. - f[7]/36. + f[8]/36. + f[9]/24. - f[10]/12. + f[11]/72. - f[16]/12. - f[17]/8. + f[18]/4. - f[19]/24. + f[20]/6. + f[21]/4. - f[22]/2. + f[23]/12. - f[24]/12. - f[25]/8. + f[26]/4. - f[27]/24. + f[32]/12. + f[33]/8. - f[34]/4. + f[35]/24. - f[36]/6. - f[37]/4. + f[38]/2. - f[39]/12. + f[40]/12. + f[41]/8. - f[42]/4. + f[43]/24. - f[48]/36. - f[49]/24. + f[50]/12. - f[51]/72. + f[52]/18. + f[53]/12. - f[54]/6. + f[55]/36. - f[56]/36. - f[57]/24. + f[58]/12. - f[59]/72.;
  a[58] = -f[0]/24. + f[1]/12. - f[2]/24. + f[4]/12. - f[5]/6. + f[6]/12. - f[8]/24. + f[9]/12. - f[10]/24. + f[16]/8. - f[17]/4. + f[18]/8. - f[20]/4. + f[21]/2. - f[22]/4. + f[24]/8. - f[25]/4. + f[26]/8. - f[32]/8. + f[33]/4. - f[34]/8. + f[36]/4. - f[37]/2. + f[38]/4. - f[40]/8. + f[41]/4. - f[42]/8. + f[48]/24. - f[49]/12. + f[50]/24. - f[52]/12. + f[53]/6. - f[54]/12. + f[56]/24. - f[57]/12. + f[58]/24.;
  a[59] = f[0]/72. - f[1]/24. + f[2]/24. - f[3]/72. - f[4]/36. + f[5]/12. - f[6]/12. + f[7]/36. + f[8]/72. - f[9]/24. + f[10]/24. - f[11]/72. - f[16]/24. + f[17]/8. - f[18]/8. + f[19]/24. + f[20]/12. - f[21]/4. + f[22]/4. - f[23]/12. - f[24]/24. + f[25]/8. - f[26]/8. + f[27]/24. + f[32]/24. - f[33]/8. + f[34]/8. - f[35]/24. - f[36]/12. + f[37]/4. - f[38]/4. + f[39]/12. + f[40]/24. - f[41]/8. + f[42]/8. - f[43]/24. - f[48]/72. + f[49]/24. - f[50]/24. + f[51]/72. + f[52]/36. - f[53]/12. + f[54]/12. - f[55]/36. - f[56]/72. + f[57]/24. - f[58]/24. + f[59]/72.;
  a[60] = f[1]/36. - f[5]/12. + f[9]/12. - f[13]/36. - f[17]/12. + f[21]/4. - f[25]/4. + f[29]/12. + f[33]/12. - f[37]/4. + f[41]/4. - f[45]/12. - f[49]/36. + f[53]/12. - f[57]/12. + f[61]/36.;
  a[61] = -f[0]/108. - f[1]/72. + f[2]/36. - f[3]/216. + f[4]/36. + f[5]/24. - f[6]/12. + f[7]/72. - f[8]/36. - f[9]/24. + f[10]/12. - f[11]/72. + f[12]/108. + f[13]/72. - f[14]/36. + f[15]/216. + f[16]/36. + f[17]/24. - f[18]/12. + f[19]/72. - f[20]/12. - f[21]/8. + f[22]/4. - f[23]/24. + f[24]/12. + f[25]/8. - f[26]/4. + f[27]/24. - f[28]/36. - f[29]/24. + f[30]/12. - f[31]/72. - f[32]/36. - f[33]/24. + f[34]/12. - f[35]/72. + f[36]/12. + f[37]/8. - f[38]/4. + f[39]/24. - f[40]/12. - f[41]/8. + f[42]/4. - f[43]/24. + f[44]/36. + f[45]/24. - f[46]/12. + f[47]/72. + f[48]/108. + f[49]/72. - f[50]/36. + f[51]/216. - f[52]/36. - f[53]/24. + f[54]/12. - f[55]/72. + f[56]/36. + f[57]/24. - f[58]/12. + f[59]/72. - f[60]/108. - f[61]/72. + f[62]/36. - f[63]/216.;
  a[62] = f[0]/72. - f[1]/36. + f[2]/72. - f[4]/24. + f[5]/12. - f[6]/24. + f[8]/24. - f[9]/12. + f[10]/24. - f[12]/72. + f[13]/36. - f[14]/72. - f[16]/24. + f[17]/12. - f[18]/24. + f[20]/8. - f[21]/4. + f[22]/8. - f[24]/8. + f[25]/4. - f[26]/8. + f[28]/24. - f[29]/12. + f[30]/24. + f[32]/24. - f[33]/12. + f[34]/24. - f[36]/8. + f[37]/4. - f[38]/8. + f[40]/8. - f[41]/4. + f[42]/8. - f[44]/24. + f[45]/12. - f[46]/24. - f[48]/72. + f[49]/36. - f[50]/72. + f[52]/24. - f[53]/12. + f[54]/24. - f[56]/24. + f[57]/12. - f[58]/24. + f[60]/72. - f[61]/36. + f[62]/72.;
  a[63] = -f[0]/216. + f[1]/72. - f[2]/72. + f[3]/216. + f[4]/72. - f[5]/24. + f[6]/24. - f[7]/72. - f[8]/72. + f[9]/24. - f[10]/24. + f[11]/72. + f[12]/216. - f[13]/72. + f[14]/72. - f[15]/216. + f[16]/72. - f[17]/24. + f[18]/24. - f[19]/72. - f[20]/24. + f[21]/8. - f[22]/8. + f[23]/24. + f[24]/24. - f[25]/8. + f[26]/8. - f[27]/24. - f[28]/72. + f[29]/24. - f[30]/24. + f[31]/72. - f[32]/72. + f[33]/24. - f[34]/24. + f[35]/72. + f[36]/24. - f[37]/8. + f[38]/8. - f[39]/24. - f[40]/24. + f[41]/8. - f[42]/8. + f[43]/24. + f[44]/72. - f[45]/24. + f[46]/24. - f[47]/72. + f[48]/216. - f[49]/72. + f[50]/72. - f[51]/216. - f[52]/72. + f[53]/24. - f[54]/24. + f[55]/72. + f[56]/72. - f[57]/24. + f[58]/24. - f[59]/72. - f[60]/216. + f[61]/72. - f[62]/72. + f[63]/216.;
}

double HorizonStatistics::evaluate_interpolation(
  double * a, double x, double y, double z)
{
#define P2(x) x*x
#define P3(x) x*x*x

  return a[0] + a[1]*z + a[2]*P2(z) + a[3]*P3(z) + a[4]*y + a[5]*y*z
    + a[6]*y*P2(z) + a[7]*y*P3(z) + a[8]*P2(y) + a[9]*P2(y)*z + a[10]*P2(y)*P2(z)
    + a[11]*P2(y)*P3(z) + a[12]*P3(y) + a[13]*P3(y)*z + a[14]*P3(y)*P2(z)
    + a[15]*P3(y)*P3(z) + a[16]*x + a[17]*x*z + a[18]*x*P2(z) + a[19]*x*P3(z)
    + a[20]*x*y + a[21]*x*y*z + a[22]*x*y*P2(z) + a[23]*x*y*P3(z) + a[24]*x*P2(y)
    + a[25]*x*P2(y)*z + a[26]*x*P2(y)*P2(z) + a[27]*x*P2(y)*P3(z) + a[28]*x*P3(y)
    + a[29]*x*P3(y)*z + a[30]*x*P3(y)*P2(z) + a[31]*x*P3(y)*P3(z) + a[32]*P2(x)
    + a[33]*P2(x)*z + a[34]*P2(x)*P2(z) + a[35]*P2(x)*P3(z) + a[36]*P2(x)*y
    + a[37]*P2(x)*y*z + a[38]*P2(x)*y*P2(z) + a[39]*P2(x)*y*P3(z) + a[40]*P2(x)*P2(y)
    + a[41]*P2(x)*P2(y)*z + a[42]*P2(x)*P2(y)*P2(z) + a[43]*P2(x)*P2(y)*P3(z)
    + a[44]*P2(x)*P3(y) + a[45]*P2(x)*P3(y)*z + a[46]*P2(x)*P3(y)*P2(z)
    + a[47]*P2(x)*P3(y)*P3(z) + a[48]*P3(x) + a[49]*P3(x)*z + a[50]*P3(x)*P2(z)
    + a[51]*P3(x)*P3(z) + a[52]*P3(x)*y + a[53]*P3(x)*y*z + a[54]*P3(x)*y*P2(z)
    + a[55]*P3(x)*y*P3(z) + a[56]*P3(x)*P2(y) + a[57]*P3(x)*P2(y)*z
    + a[58]*P3(x)*P2(y)*P2(z) + a[59]*P3(x)*P2(y)*P3(z) + a[60]*P3(x)*P3(y)
    + a[61]*P3(x)*P3(y)*z + a[62]*P3(x)*P3(y)*P2(z) + a[63]*P3(x)*P3(y)*P3(z);
}

  
  // doing derivatives to ONLY level function:
  // r - h(\theta, \phi)
real_t HorizonStatistics::dF(int theta_i, int phi_i, int d, double x, double y, double z, double r)
{
  double res;
  if(theta_i + 1 >= 2 * n_theta || theta_i < 0)
    TBOX_ERROR("Theta_i is out of the range!\b");
  
  double dphi = 2.0 * PI / (double) n_phi / 2.0;
  double dtheta = PI / (double) n_theta / 2.0;
  
  double dtheta_dh = (ah_radius[theta_i + 1][phi_i] - ah_radius[theta_i][phi_i]) / (dtheta);
  double dphi_dh = (ah_radius[theta_i][(phi_i+1)%(2*n_phi)]
                    - ah_radius[theta_i][(phi_i)%(2*n_phi)]) / (dphi);
  if(d == 1)
    res = x/r - (x*z/(pw3(r)*sqrt(1-pw2(z/r))))*dtheta_dh + (y/(pw2(x)+pw2(y)))*dphi_dh;
  else if(d == 2)
    res = y/r - (y*z/(pw3(r)*sqrt(1-pw2(z/r))))*dtheta_dh - (x/(pw2(x)+pw2(y)))*dphi_dh;
  else if(d == 3)
    res = z/r + sqrt(pw2(x)+pw2(y)) / pw2(r) * dtheta_dh;
  else
    TBOX_ERROR("Direction "<<d<< "is not recognized\n");
  return res;
}
  
real_t HorizonStatistics::findRadius(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double theta_0, double phi_0)
{
  return getRadius(theta_0, phi_0);
}
void HorizonStatistics::set_G_values(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double theta, double phi, int theta_i, int phi_i, double r, KillingData *kd, BSSN * bssn)
{

}


void HorizonStatistics::set_kd_values(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double theta, double phi, int theta_i, int phi_i, double r, KillingData *kd, BSSN * bssn)
{

}

real_t HorizonStatistics::ev_k_theta_dphi(KillingData *kd, int theta_i, int phi_i)
{
  return kd->Gs112 * k_theta[theta_i][phi_i] + kd->Gs212 * k_phi[theta_i][phi_i]
    - k_L[theta_i][phi_i] * sqrt(kd->q11 * kd->q22 - kd->q12 * kd->q12);
}

real_t HorizonStatistics::ev_k_phi_dphi(KillingData *kd, int theta_i, int phi_i)
{
  return kd->Gs122 * k_theta[theta_i][phi_i] + kd->Gs222 * k_phi[theta_i][phi_i];
}

real_t HorizonStatistics::ev_k_L_dphi(KillingData *kd, int theta_i, int phi_i)
{
  return 0.5 * kd->R * sqrt(kd->q11 * kd->q22 - kd->q12 * kd->q12)
    * (kd->qi11 * k_theta[theta_i][phi_i] + kd->qi12 * k_phi[theta_i][phi_i]);
}

real_t HorizonStatistics::ev_k_theta_dtheta(KillingData *kd, int theta_i, int phi_i)
{
  return kd->Gs111 * k_theta[theta_i][phi_i] + kd->Gs211 * k_phi[theta_i][phi_i];
}

real_t HorizonStatistics::ev_k_phi_dtheta(KillingData *kd, int theta_i, int phi_i)
{
  return kd->Gs112 * k_theta[theta_i][phi_i] + kd->Gs212 * k_phi[theta_i][phi_i]
    + k_L[theta_i][phi_i] * sqrt(kd->q11 * kd->q22 - kd->q12 * kd->q12);
}

real_t HorizonStatistics::ev_k_L_dtheta(KillingData *kd, int theta_i, int phi_i)
{
  return -0.5 * kd->R * sqrt(kd->q11 * kd->q22 - kd->q12 * kd->q12)
        * (kd->qi12 * k_theta[theta_i][phi_i] + kd->qi22 * k_phi[theta_i][phi_i]);
    //* (kd->qi11 * k_theta[theta_i][phi_i] + kd->qi12 * k_phi[theta_i][phi_i]);
}

real_t HorizonStatistics::getRadius(double theta_i, double phi_i)
{
  real_t x = 1.0 * cos(phi_i) * sin(theta_i) + origin[0];
  real_t y = 1.0 * sin(phi_i) * sin(theta_i) + origin[1];
  real_t z = 1.0 * cos(theta_i) + origin[2];
  real_t res = 0;
  horizon->AHFinderDirect_radius_in_direction(
    horizon_id, 1, &x, &y, &z, &res);

  return res;

}
  
void HorizonStatistics::transportKillingTheta(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t phi_i, double k_theta_0, double k_phi_0, double k_L_0, BSSN * bssn)
{
}

void HorizonStatistics::transportKillingPhi(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t theta_i, idx_t phi_f, double k_theta_0, double k_phi_0, double k_L_0, BSSN * bssn)
{
}

void HorizonStatistics::initG(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn)
{

}



void HorizonStatistics::initGridding(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{

}

void HorizonStatistics::findM(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, double x[], BSSN *bssn)
{
  //  using namespace std::literals::complex_literals;

  transportKillingPhi(hierarchy, n_theta/2, n_phi, 1, 0, 0, bssn);

  Eigen::Matrix3d M;
  
  M(0, 0) = k_theta[n_theta/2][0];
  M(1, 0) = k_phi[n_theta/2][0];
  M(2, 0) = k_L[n_theta/2][0];


  transportKillingPhi(hierarchy, n_theta/2, n_phi, 0, 1, 0, bssn);

  M(0, 1) = k_theta[n_theta/2][0];
  M(1, 1) = k_phi[n_theta/2][0];
  M(2, 1) = k_L[n_theta/2][0];


  
  transportKillingPhi(hierarchy, n_theta/2, n_phi, 0, 0, 1, bssn);

  M(0, 2) = k_theta[n_theta/2][0]; 
  M(1, 2) = k_phi[n_theta/2][0];
  M(2, 2) = k_L[n_theta/2][0];



   tbox::pout<<"\n";
  tbox::pout << "Here is the matrix m:\n" << M << "\n";

  // solve the eigenvalue equation to get 3 eigenvalues;

  Eigen::EigenSolver<Eigen::Matrix3d> eigensolver(M);

  if (eigensolver.info() != Eigen::Success)
    TBOX_ERROR("Matrix does not have solution with identity eigenvalue!\n");

  Eigen::EigenSolver< Eigen::Matrix3d >::EigenvalueType e_val = eigensolver.eigenvalues();
  Eigen::EigenSolver< Eigen::Matrix3d >::EigenvectorsType e_vec = eigensolver.eigenvectors();

   tbox::pout<<"\n Eigenvalues are "<<e_val<<"\n";
  
  double dis_to_I = INF;
  int identity_idx=-1;
  for(int i = 0; i < 3; i++)
    if(pw2((e_val(i).real() - 1.0)) + pw2(e_val(i).imag()) < dis_to_I)
    {
      dis_to_I = pw2(e_val(i).real() - 1.0) + pw2(e_val(i).imag());
      identity_idx = i;
    }

  x[0] = e_vec(0, identity_idx).real();
  x[1] = e_vec(1, identity_idx).real();
  x[2] = e_vec(2, identity_idx).real();
  

   tbox::pout<<"Solution with identity eigenvalue is ("
         << x[0]<<" "<<x[1]<<" "<<x[2]<<")\n";
  
  
}

void HorizonStatistics::set_norm_values(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double theta, double phi, int theta_i, int phi_i, double r, KillingData *kd, BSSN * bssn)
{
}

real_t HorizonStatistics::interp_k_theta(
  double theta, double phi)
{
  if(theta > PI) theta = 2* PI - theta;
  if(theta < 0) theta = -theta;

  if(phi > 2.0 * PI) phi -= 2.0 * PI;
  if(phi < 0) phi = 2.0 * PI + phi;

  double d_theta = PI / (double)n_theta;
  double d_phi = 2.0 * PI / (double)n_phi;
  int theta_i0 = floor((theta) / d_theta - 0.5);
  int phi_i0 = floor((phi) / d_phi - 0.5);

  phi_i0 = (phi_i0 + n_phi)%n_phi;

  if(phi_i0 >= n_phi || phi_i0 < 0 || theta_i0 >= n_theta || theta_i0 < 0)
  {
    TBOX_ERROR("Step is too large!\n "<<theta<<" "<<phi<<" "<<theta_i0<<" "<<phi_i0<<"\n");
  }

  double res = 0;
  for(int theta_i = theta_i0; theta_i != (theta_i0+2)%n_theta; theta_i = (theta_i +1 )%n_theta)
  {
    double theta0 = PI * ((double) theta_i +0.5) / (double)n_theta;

    double temp = 0;
    for(int phi_i = phi_i0; phi_i != (phi_i0+2)%n_phi; phi_i = (phi_i +1 )%n_phi)
    {
      double phi0 = 2.0 * PI * ((double) phi_i +0.5) / (double)n_phi;

      double phi_dist = phi0 - phi;
      if(phi_dist > d_phi + EPS)
        phi_dist -= 2.0 * PI;
      else if(phi_dist < -d_phi - EPS)
        phi_dist += 2.0 * PI;
      if(phi_dist > d_phi + EPS || phi_dist < -d_phi - EPS)
        TBOX_ERROR("Error in calculating phi_dist "<<phi0<<" "<<phi<<" "<<phi_dist<<" "<<d_phi<<"\n");
      temp += k_theta[theta_i][phi_i] *        
        ( (d_phi - (phi_dist ) * (2.0 * ( (phi_i - phi_i0 + n_phi)%n_phi ) -1.0)) / d_phi);
    }
    res += temp *
      ( (d_theta - (theta0 - theta ) * (2.0 * (theta_i - theta_i0) -1.0)) / d_theta);
  }
  
  return res;
}

real_t HorizonStatistics::interp_k_phi(
  double theta, double phi)
{
  if(theta > PI) theta = 2* PI - theta;
  if(theta < 0) theta = -theta;

  if(phi > 2.0 * PI) phi -= 2.0 * PI;
  if(phi < 0) phi = 2.0 * PI + phi;

  double d_theta = PI / (double)n_theta;
  double d_phi = 2.0 * PI / (double)n_phi;
  int theta_i0 = floor((theta) / d_theta - 0.5);
  int phi_i0 = floor((phi) / d_phi - 0.5);

  phi_i0 = (phi_i0 + n_phi)%n_phi;
  
  if(phi_i0 > n_phi || phi_i0 < 0 || theta_i0 > n_theta || theta_i0 < 0)
  {
    TBOX_ERROR("EEEE "<<theta<<" "<<phi<<" "<<theta_i0<<" "<<phi_i0<<"\n");
  }
  double res = 0;

  for(int theta_i = theta_i0; theta_i != (theta_i0+2)%n_theta;  theta_i = (theta_i +1 )%n_theta)
  {
    double theta0 = PI * ((double) theta_i +0.5) / (double)n_theta;

    double temp = 0;
    for(int phi_i = phi_i0; phi_i != (phi_i0+2)%n_phi; phi_i = (phi_i +1 )%n_phi)
    {
      double phi0 = 2.0 * PI * ((double) phi_i +0.5) / (double)n_phi;

      double phi_dist = phi0 - phi;
      if(phi_dist > d_phi + EPS)
        phi_dist -= 2.0 * PI;
      else if(phi_dist < -d_phi - EPS)
        phi_dist += 2.0 * PI;
      if(phi_dist > d_phi + EPS || phi_dist < -d_phi - EPS)
        TBOX_ERROR("Error in calculating phi_dist "<<phi0<<" "<<phi<<" "<<phi_dist<<" "<<d_phi<<"\n");
      temp += k_phi[theta_i][phi_i] *        
        ( (d_phi - (phi_dist ) * (2.0 * ( (phi_i - phi_i0 + n_phi)%n_phi ) -1.0)) / d_phi);

    }
    res += temp *
      ( (d_theta - (theta0 - theta ) * (2.0 * (theta_i - theta_i0) -1.0)) / d_theta);
  }
  
  return res;

}

void HorizonStatistics::normKilling()
{
  double c = getNormFactor();

   tbox::pout<<"Norm factor "<<c<<"\n";
  for(int i = 0 ; i < n_theta; i++)
    for(int j = 0; j < n_phi; j ++)
    {
      k_theta[i][j] *= c;
      k_phi[i][j] *= c;
    }
}

real_t HorizonStatistics::getNormFactor()
{
  double theta = PI/2, phi = 0;

  real_t dtheta = PI / (double)n_theta;
  real_t pre_phi = 0;

  real_t max_abs = 0;

  for(int i = 0; i < n_theta; i ++)
    for(int j = 0; j < n_phi; j ++)
      max_abs = std::max(max_abs, std::max(fabs(k_phi[i][j]), fabs(k_theta[i][j])));
  
  real_t dt = 0.01 * dtheta / max_abs;

  tbox::pout<<"Time interval for process of getting norm factor is "
            <<dt<<"\n";

  
  real_t t = 0;
  int cnt = 0;
  // advance theta and phi
  while( (SIGN(pre_phi) * SIGN(phi) >= 0)
         && (SIGN(2.0 * PI - pre_phi) * SIGN(2.0 * PI - phi) >= 0 )
         &&(SIGN(-2.0 * PI - pre_phi) * SIGN(-2.0 * PI - phi) >= 0 )) 
  {
    pre_phi = phi;

    double theta_0 = theta;
    double phi_0 = phi;
    
    real_t k1_theta = interp_k_theta(theta, phi);
    real_t k1_phi = interp_k_phi(theta, phi);

    theta = theta_0 + k1_theta * dt / 2.0, phi = phi_0 + k1_phi * dt / 2.0;
    real_t k2_theta = interp_k_theta(theta , phi);
    real_t k2_phi = interp_k_phi(theta , phi);

    theta = theta_0 + k2_theta * dt / 2.0, phi = phi_0 + k2_phi * dt / 2.0;
    real_t k3_theta = interp_k_theta(theta, phi);
    real_t k3_phi = interp_k_phi(theta, phi);

    theta = theta_0 + k3_theta * dt, phi = phi_0 + k3_phi * dt;
    real_t k4_theta = interp_k_theta(theta , phi);
    real_t k4_phi = interp_k_phi(theta , phi);

    theta = theta_0 + dt * (k1_theta + 2.0 * k2_theta + 2.0 * k3_theta + k4_theta) / 6.0;
    phi = phi_0 + dt * (k1_phi + 2.0 * k2_phi + 2.0 * k3_phi + k4_phi) / 6.0;
    t += dt;
    cnt++;
    if(cnt >= 1000000)
    {
      tbox::pout<<"Warning!!! Cannot find suitable normalize factor\n";
      break;
    }
  }

  //  tbox::pout<<"When \phi gets back, the corresponding theta is "<<theta<<"\n";
  return t/ (2.0 * PI) ;
}

// here must be killing "vector" (not one form)
real_t HorizonStatistics::angularMomentum(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn)
{
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  real_t dtheta = PI / (double)n_theta;
  real_t dphi = 2.0 * PI / (double)n_phi;
  real_t res = 0;
  for(int theta_i = 0; theta_i < n_theta; theta_i++)
  {
    for(int phi_i = 0; phi_i < n_phi; phi_i++)
    {
      //      std::cout<<theta_i<<" "<<phi_i<<" "<<res<<"\n";
      double theta = PI * ((double) theta_i +0.5) / (double)n_theta;
      double phi = 2.0 * PI * ((double) phi_i +0.5) / (double)n_phi;

      real_t st = sin(theta);
      real_t ct = cos(theta);
      real_t sp = sin(phi);
      real_t cp = cos(phi);
      
      real_t r = getRadius(theta, phi);
      
      KillingData kd = {0};

      set_norm_values(hierarchy, theta, phi, theta_i, phi_i,
                      getRadius(theta, phi), &kd, bssn);
      

      real_t k1 = r * (k_theta[theta_i][phi_i] * ct * cp
                                                   - k_phi[theta_i][phi_i] * sp * st);
      real_t k2 = r * (k_theta[theta_i][phi_i] * ct * sp
                                                   + k_phi[theta_i][phi_i] * cp * st);
      real_t k3 = - r * k_theta[theta_i][phi_i] * st;

      
      real_t s1 = (kd.mi11 * kd.d1F + kd.mi12 * kd.d2F + kd.mi13 * kd.d3F)
        / (sqrt((kd.mi11 * kd.d1F * kd.d1F + kd.mi22 * kd.d2F * kd.d2F + kd.mi33 * kd.d3F *kd.d3F
                 + 2.0 * (kd.mi12 * kd.d1F * kd.d2F + kd.mi13 * kd.d1F * kd.d3F + kd.mi23 * kd.d2F * kd.d3F))));
      real_t s2 = (kd.mi21 * kd.d1F + kd.mi22 * kd.d2F + kd.mi23 * kd.d3F)
        / (sqrt((kd.mi11 * kd.d1F * kd.d1F + kd.mi22 * kd.d2F * kd.d2F + kd.mi33 * kd.d3F *kd.d3F
                 + 2.0 * (kd.mi12 * kd.d1F * kd.d2F + kd.mi13 * kd.d1F * kd.d3F + kd.mi23 * kd.d2F * kd.d3F))));
      real_t s3 = (kd.mi31 * kd.d1F + kd.mi32 * kd.d2F + kd.mi33 * kd.d3F)
        / (sqrt((kd.mi11 * kd.d1F * kd.d1F + kd.mi22 * kd.d2F * kd.d2F + kd.mi33 * kd.d3F *kd.d3F
                 + 2.0 * (kd.mi12 * kd.d1F * kd.d2F + kd.mi13 * kd.d1F * kd.d3F + kd.mi23 * kd.d2F * kd.d3F))));

      
      real_t K11 = kd.K11, K12 = kd.K12, K13 = kd.K13;
      real_t K22 = kd.K22, K23 = kd.K23, K33 = kd.K33;

      
      kd = {0};

      set_kd_values(hierarchy, theta, phi, theta_i*2, phi_i*2,
                    r, &kd, bssn);

      double det = kd.q11 * kd.q22 - kd.q12 * kd.q12;

      res += (k1 * s1 * K11 + k2 * s2 * K22 + k3 * s3 * K33
        + k1 * s2 * K12 + k1 * s3 * K13 + k2 * s3 * K23
              + k2 * s1 * K12 + k3 * s1 * K13 + k3 * s2 * K23) * sqrt(det) * dtheta * dphi;

      
      mpi.Barrier();
      if (mpi.getSize() > 1 ) {
        mpi.Bcast(&res, 1, MPI_DOUBLE, cur_mpi_rank);
      }
      mpi.Barrier();
    }
  }

  return res / 8.0 / PI;
}

void HorizonStatistics::convertToVector(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn)
{
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  for(int theta_i = 0; theta_i < n_theta; theta_i++)
  {
    for(int phi_i = 0; phi_i < n_phi; phi_i++)
    {
      double k_theta0 = k_theta[theta_i][phi_i];
      double k_phi0 = k_phi[theta_i][phi_i];
      double theta = PI * ((double) theta_i +0.5) / (double)n_theta;
      double phi = 2.0 * PI * ((double) phi_i +0.5) / (double)n_phi;

      KillingData kd = {0};

      set_kd_values(hierarchy, theta, phi, theta, phi, getRadius(theta, phi), &kd, bssn);
      
      k_theta[theta_i][phi_i] = kd.qi11 * k_theta0 + kd.qi12 * k_phi0;
      k_phi[theta_i][phi_i] = kd.qi12 * k_theta0 + kd.qi22 * k_phi0;
      
      mpi.Barrier();
      if (mpi.getSize() > 1 ) {
        mpi.Bcast(&k_theta[theta_i][phi_i], 1, MPI_DOUBLE, cur_mpi_rank);
        mpi.Bcast(&k_phi[theta_i][phi_i], 1, MPI_DOUBLE, cur_mpi_rank);
      }
      mpi.Barrier();

    }
  }


}

real_t HorizonStatistics::area(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn)
{
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());

  real_t dtheta = PI / (double)n_theta;
  real_t dphi = 2.0 * PI / (double)n_phi;
  real_t res = 0;
  for(int theta_i = 0; theta_i < n_theta; theta_i++)
  {
    for(int phi_i = 0; phi_i < n_phi; phi_i++)
    {
      double theta = PI * ((double) theta_i +0.5) / (double)n_theta;
      double phi = 2.0 * PI * ((double) phi_i +0.5) / (double)n_phi;

      real_t r = getRadius(theta, phi);
      
      KillingData kd = {0};

      set_kd_values(hierarchy, theta, phi, theta_i*2, phi_i*2,
                    r, &kd, bssn);

      double det = kd.q11 * kd.q22 - kd.q12 * kd.q12;

      res += sqrt(det) * dtheta * dphi;

      mpi.Barrier();
      if (mpi.getSize() > 1 ) {
        mpi.Bcast(&res, 1, MPI_DOUBLE, cur_mpi_rank);
      }
      mpi.Barrier();
    }
  }
  return res;

}


void HorizonStatistics::findKilling(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, BSSN * bssn, int horizon_id_in, int step)
{
  if(step % horizon->find_every != 0)
    return;
  horizon_id = horizon_id_in;

  horizon->AHFinderDirect_local_coordinate_origin(
    horizon_id, &origin[0],&origin[1],&origin[2]);
  
  tbox::pout<<"Starting the process of finding Killing vectors for horions: "<<horizon_id<<"\n!";

  initGridding(hierarchy);

  initG(hierarchy, bssn);

  double angular_m = 0;

  
  if(non_zero_angular_momentum)
  {
    double x[3] = {0};
    // Finding transport matrix and initializing eigen vector
    findM(hierarchy, x, bssn);
  
    transportKillingPhi(hierarchy, n_theta/2, n_phi - 1, x[0], x[1], x[2], bssn);

    for(int i = 0; i < n_phi; i++)
    {
      transportKillingTheta(
        hierarchy, i, k_theta[n_theta/2][i], k_phi[n_theta/2][i], k_L[n_theta/2][i], bssn);
    }

    // transportKillingTheta(hierarchy, 0, x[0], x[1], x[2], bssn);

    // for(int i = 0; i < n_theta; i++)
    // {
    //   transportKillingPhi(
    //     hierarchy, i, n_phi-1, k_theta[i][0], k_phi[i][0], k_L[i][0], bssn);
    // }
  
    convertToVector(hierarchy, bssn);

    normKilling();
    
    angular_m = angularMomentum(hierarchy, bssn);
    tbox::pout<<"Angular momentum is "<<angular_m<<"\n";
  }
  
  double a = area(hierarchy, bssn);

  double R_Delta = sqrt(a / (4.0 * PI));

  double bare_mass = sqrt(a / (4.0 * PI)) / 2.0;
  double mass = sqrt(pw2(pw2(R_Delta)) + 4.0 * pw2(angular_m)) / (2.0 * R_Delta);

  tbox::pout<<"Mass is "<<mass<<" irreducible mass (areal) is "<<bare_mass<<"\n";
  
}
  
}
