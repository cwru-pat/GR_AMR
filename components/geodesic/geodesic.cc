#include "geodesic.h"
#include "../../cosmo_includes.h"
#include "../../utils/math.h"


using namespace SAMRAI;

namespace cosmo
{
Geodesic::Geodesic(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  const tbox::Dimension& dim_in,
  std::shared_ptr<tbox::Database> database_in,
  std::ostream* l_stream_in,
  real_t KO_damping_coefficient_in,
  int weight_idx_in):
  cosmo_geodesic_db(database_in),
  lstream(l_stream_in),
  dim(dim_in),
  weight_idx(weight_idx_in),
  ghost_width(cosmo_geodesic_db->getIntegerWithDefault("ghost_width", 2)),
  p0(cosmo_geodesic_db->getDoubleWithDefault("p0", 1)),
  save_metric(cosmo_geodesic_db->getBoolWithDefault("save_metric", false)),
  pc(new pdat::IndexVariable<ParticleContainer, pdat::CellGeometry>(
      dim, "particle"))
{
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

  std::shared_ptr<hier::VariableContext> context_scratch(
    variable_db->getContext("SCRATCH"));
  std::shared_ptr<hier::VariableContext> context_active(
    variable_db->getContext("ACTIVE"));
  std::shared_ptr<hier::VariableContext> context_previous(
    variable_db->getContext("PREVIOUS"));
  std::shared_ptr<hier::VariableContext> context_k1(
    variable_db->getContext("RK_K1"));
  std::shared_ptr<hier::VariableContext> context_k2(
    variable_db->getContext("RK_K2"));
  std::shared_ptr<hier::VariableContext> context_k3(
    variable_db->getContext("RK_K3"));
  std::shared_ptr<hier::VariableContext> context_k4(
    variable_db->getContext("RK_K4"));

  // new context as buffer storing particles that going to
  // coearser level
  std::shared_ptr<hier::VariableContext> context_down_stream_buffer(
    variable_db->getContext("down_stream_buffer"));

    
  pc_idx = variable_db->registerVariableAndContext(
    pc, context_active, hier::IntVector(dim, ghost_width));

  pc_s_idx = variable_db->registerVariableAndContext(
    pc, context_scratch, hier::IntVector(dim, ghost_width));


  pc_d_buffer_idx = variable_db->registerVariableAndContext(
    pc, context_down_stream_buffer, hier::IntVector(dim, ghost_width));


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

  for(int i = 0 ; i < 3; i++)
    L[i] = domain_upper[i] - domain_lower[i];


  if(save_metric && PARTICLE_REAL_PROPERTIES < 7)
    TBOX_ERROR("Need store metric into particle state but no enough space!");

  if(ghost_width + 3 > GHOST_WIDTH)
    TBOX_ERROR("There is NO enough ghost with for fields to do splines interpolation for particles!");
}

void Geodesic::compute_tricubic_coeffs(double *a, double *f)
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

double Geodesic::evaluate_interpolation(
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

  
// initializing all particles
void Geodesic::initAll(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  
  std::string init_type = cosmo_geodesic_db->getStringWithDefault("init_type","");
        
  if(init_type == "Schwarzchild")
  {
    geodesic_ic_Schwarzchild_test(
      hierarchy, cosmo_geodesic_db);
  }
  else if(init_type == "2rays")
  {
    geodesic_ic_face_null_test(
      hierarchy, cosmo_geodesic_db);
  }
  else
    TBOX_ERROR("Unsupported null geodesic initial type!");
}
  
void Geodesic::initPData(
  const std::shared_ptr<hier::Patch> & patch)
{
  pc_pdata =
    SAMRAI_SHARED_PTR_CAST<pdat::IndexData<ParticleContainer, pdat::CellGeometry>,
               hier::PatchData>(
                 patch->getPatchData(pc_idx));

  pc_s_pdata =
    SAMRAI_SHARED_PTR_CAST<pdat::IndexData<ParticleContainer, pdat::CellGeometry>,
               hier::PatchData>(
                 patch->getPatchData(pc_s_idx));

  pc_d_buffer_pdata =
    SAMRAI_SHARED_PTR_CAST<pdat::IndexData<ParticleContainer, pdat::CellGeometry>,
               hier::PatchData>(
                 patch->getPatchData(pc_d_buffer_idx));
  
}

void Geodesic::allocParticles(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int ln)
{
  std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
  level->allocatePatchData(pc_idx);
  level->allocatePatchData(pc_s_idx);
  level->allocatePatchData(pc_d_buffer_idx);
}

void Geodesic::insertPatchParticles(
  const std::shared_ptr<hier::Patch> & patch, int src_id, int dst_id)
{
  
  std::shared_ptr<pdat::IndexData<ParticleContainer,
                                    pdat::CellGeometry> > src_pdata(
                                      SAMRAI_SHARED_PTR_CAST<pdat::IndexData<ParticleContainer, pdat::CellGeometry>,
                                      hier::PatchData>(
                                        patch->getPatchData(src_id)));
  std::shared_ptr<pdat::IndexData<ParticleContainer,
                                    pdat::CellGeometry> > dst_pdata(
                                      SAMRAI_SHARED_PTR_CAST<pdat::IndexData<ParticleContainer, pdat::CellGeometry>,
                                      hier::PatchData>(
                                        patch->getPatchData(dst_id)));

  const hier::Box box = src_pdata->getGhostBox();

  pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iter(*src_pdata, true);
  pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iterend(*src_pdata, false);

  for(; iter != iterend; iter++)
  {
    // get pointer to this particle at Index
    ParticleContainer & src_p = *iter;
    hier::Index ic (iter.getIndex());
    ParticleContainer *dst_p = dst_pdata->getItem(ic);

    if(dst_p == NULL)
    {
      ParticleContainer *p = new ParticleContainer(ic);
      dst_pdata->appendItem(ic, *p);
      dst_p = dst_pdata->getItem(ic);
    }
    dst_p->p_list.insert(dst_p->p_list.end(), src_p.p_list.begin(), src_p.p_list.end());
    
  }

  
}

void Geodesic::clearParticlesLivingInGhostCells(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
    clearParticlesLivingInGhostCells(hierarchy, ln);
}


void Geodesic::clearParticlesLivingInGhostCells(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int ln)
{
  if(ln > hierarchy->getNumberOfLevels() - 1)
    return;
  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  // loop over patches on level
  for (hier::PatchLevel::iterator ip(level->begin());
       ip != level->end(); ++ip)
  {
    std::shared_ptr<hier::Patch> patch(*ip);

    initPData(patch);
    pc_pdata->removeGhostItems();
    pc_s_pdata->removeGhostItems();
  }
  
}


void Geodesic::insertLevelParticles(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int ln, int src_id, int dst_id)
{
  std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end(); ++ip)
    {
      std::shared_ptr<hier::Patch> patch(*ip);
      insertPatchParticles(patch, src_id, dst_id);
    }
}

// clear particles for all hierarchy
void Geodesic::clearParticles(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int idx)
{
  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
  {
    clearParticles(hierarchy, ln, idx);
  }

}

// clear particles for level ln
void Geodesic::clearParticles(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int ln, int idx)
{
  
  std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
    //              std::cout<<ln<<" "<<patch->getBox()<<"\n"<<std::flush;

    clearParticles(patch, idx);
  }  

}

// clear particles for a patch
void Geodesic::clearParticles(
  const std::shared_ptr<hier::Patch> & patch, int idx)
{
  std::shared_ptr<pdat::IndexData<ParticleContainer,
                                  pdat::CellGeometry> > pdata(
                                    SAMRAI_SHARED_PTR_CAST<pdat::IndexData<ParticleContainer, pdat::CellGeometry>,
                                    hier::PatchData>(
                                      patch->getPatchData(idx)));
  
  // delete all items in this patch data
  pdata->removeAllItems();
}


void Geodesic::preAdvance(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int ln)
{
  std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
  
  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  
  TBOX_ASSERT(grid_geometry_);
  
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;


  // before advancing, should fill ghost cells at the same level
  // note these particles are REAL particles living in ghost cells
  // so the filling should be from real container to real container
 {  
   std::shared_ptr<hier::RefineOperator> refine_op
     = grid_geometry.
     lookupRefineOperator(pc, "PARTICLE_REFINE");
  

   xfer::RefineAlgorithm refiner;
  
   refiner.registerRefine(pc_idx,
                          pc_idx,                        
                          pc_idx,                        
                          refine_op);

   std::shared_ptr<xfer::RefineSchedule> refine_schedule;
  
   refine_schedule =
     refiner.createSchedule(
       level,
       NULL);

   refine_schedule->fillData(0.0);

 }
 // std::cout<<"After pre advancing\n";
 // printAll(hierarchy, pc_idx);
}
  
void Geodesic::RKEvolvePatch(
  const std::shared_ptr<hier::Patch> & patch, BSSN *bssn, real_t dt)
{
  initPData(patch);
  bssn->initPData(patch);
  bssn->initMDA(patch);

  const hier::Box& box = patch->getBox();
  
  const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(  
    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  //initialize dx for each patch
  const real_t * dx = &(patch_geom->getDx())[0];

  hier::Box ghost_box = pc_pdata->getGhostBox();
  // do not need ghost box to be that wide
  hier::Box effective_ghost_box = ghost_box;
  if(ghost_width > 1)
    effective_ghost_box.grow(hier::Index (1-ghost_width, 1-ghost_width, 1-ghost_width));
  hier::Box fields_ghost_box = bssn->DIFFchi_a_pdata->getGhostBox();
  pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iter(*pc_pdata, true);
  pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iterend(*pc_pdata, false);

  double ghost_box_phys_lower[3], ghost_box_phys_upper[3];

  for(int i = 0; i < 3; i ++)
  {
    ghost_box_phys_lower[i] = (effective_ghost_box.lower()[i] * dx[i]);
    ghost_box_phys_upper[i] = ( (effective_ghost_box.upper()[i] + 1.0) * dx[i]);
  }
    
  // pdat::CellIterator icend(pdat::CellGeometry::end(ghost_box));
  // for (pdat::CellIterator ic(pdat::CellGeometry::begin(ghost_box));
  //      ic != icend; ++ic)
  for(;iter != iterend; iter++)
  {
    //    hier::Index idx(*ic);
    ParticleContainer & id = *iter;

    double shift[3] = {0};
    

    for(std::list<RKParticle>::iterator it=id.p_list.begin();
        it != id.p_list.end(); it++)
    {
      GeodesicData gd = {0};

      for(int i = 0 ; i < 3; i ++)
      {
        if(ghost_box_phys_lower[i] > ((*it).x_a[i])) 
          shift[i] = ceil( (double)(ghost_box_phys_lower[i] - (*it).x_a[i]) / round(L[i]));
        if(ghost_box_phys_upper[i] < ((*it).x_a[i])) 
          shift[i] = -ceil( (double)(-ghost_box_phys_upper[i] + (*it).x_a[i]) / round(L[i]));

      }
      int i0 = floor(((*it).x_a[0] + shift[0] * L[0] - domain_lower[0] ) / dx[0] );
      int j0 = floor(((*it).x_a[1] + shift[1] * L[1] - domain_lower[1] ) / dx[1] );
      int k0 = floor(((*it).x_a[2] + shift[2] * L[2] - domain_lower[2] ) / dx[2] );

      // if particle is outside the ghost box anytime
      // during the RK advance, do not advance it
      hier::Index temp_idx(i0, j0, k0);

      if(!effective_ghost_box.contains(temp_idx))
      {
        for(int i = 0; i < PARTICLE_NUMBER_OF_STATES; i++)
          (*it).x_c[i] = 0;
        continue;
      }
      
      for(int i = 0; i < 8; i++)
        for(int j = 0; j < 8; j++)
          for(int k = 0; k < 8; k++)
          {
            int i0 = floor(((*it).x_a[0] + shift[0] * L[0] - domain_lower[0] ) / dx[0] -0.5);
            int j0 = floor(((*it).x_a[1] + shift[1] * L[1] - domain_lower[1] ) / dx[1] -0.5);
            int k0 = floor(((*it).x_a[2] + shift[2] * L[2] - domain_lower[2] ) / dx[2] -0.5);

            hier::Index temp_idx(i0-3+i, j0-3+j, k0-3+k);
            
            if(!fields_ghost_box.contains(temp_idx))
              TBOX_ERROR("Particle at "<<hier::Index(i0, j0, k0)
                         <<"is not contained in the ghostbox "<<fields_ghost_box<<"\n");
          }
      set_gd_values(patch, (*it).x_a, &gd, bssn, dx, shift);
      RKEvolveParticle((*it), gd, dt);

      if(PARTICLE_REAL_PROPERTIES > 0)
      {
        (*it).rp[0] = gd.p0;
        if(save_metric && PARTICLE_REAL_PROPERTIES > 6)
        {
          (*it).rp[1] = gd.m11;
          (*it).rp[2] = gd.m12;
          (*it).rp[3] = gd.m13;
          (*it).rp[4] = gd.m22;
          (*it).rp[5] = gd.m23;
          (*it).rp[6] = gd.m33;
        }
      }
    }
    
  }

}

void Geodesic::RKEvolveParticle(
  RKParticle &p, GeodesicData &gd, double dt)
{
  p.x_c[0] = dt *
    (-gd.beta1 +
     (gd.mi11 * gd.q1 + gd.mi12 * gd.q2 + gd.mi13 * gd.q3) / gd.p0);
  p.x_c[1] = dt *
    (-gd.beta2 +
     (gd.mi12 * gd.q1 + gd.mi22 * gd.q2 + gd.mi23 * gd.q3) / gd.p0);
  p.x_c[2] = dt *
    (-gd.beta3 +
     (gd.mi13 * gd.q1 + gd.mi23 * gd.q2 + gd.mi33 * gd.q3) / gd.p0);

  p.x_c[3] = dt * (
    -gd.p0 * gd.alpha * gd.d1alpha
    + (gd.q1 * gd.d1beta1 + gd.q2 * gd.d1beta2  + gd.q3 * gd.d1beta3  )
    + (gd.qi1 * gd.qi1 * gd.d1m11 + gd.qi2 * gd.qi2 * gd.d1m22 + gd.qi3 * gd.qi3 * gd.d1m33
       + 2.0*gd.qi1 * gd.qi2 * gd.d1m12 + 2.0*gd.qi1 * gd.qi3 * gd.d1m13 + 2.0*gd.qi2 * gd.qi3 * gd.d1m23)/2.0/gd.p0);

  p.x_c[4] = dt * (
    -gd.p0 * gd.alpha * gd.d2alpha
    + (gd.q1 * gd.d2beta1 + gd.q2 * gd.d2beta2  + gd.q3 * gd.d2beta3  )
    + (gd.qi1 * gd.qi1 * gd.d2m11 + gd.qi2 * gd.qi2 * gd.d2m22 + gd.qi3 * gd.qi3 * gd.d2m33
       + 2.0*gd.qi1 * gd.qi2 * gd.d2m12 + 2.0*gd.qi1 * gd.qi3 * gd.d2m13 + 2.0*gd.qi2 * gd.qi3 * gd.d2m23)/2.0/gd.p0);

  p.x_c[5] = dt * (
    -gd.p0 * gd.alpha * gd.d3alpha
    + (gd.q1 * gd.d3beta1 + gd.q2 * gd.d3beta2  + gd.q3 * gd.d3beta3  )
    + (gd.qi1 * gd.qi1 * gd.d3m11 + gd.qi2 * gd.qi2 * gd.d3m22 + gd.qi3 * gd.qi3 * gd.d3m33
       + 2.0*gd.qi1 * gd.qi2 * gd.d3m12 + 2.0*gd.qi1 * gd.qi3 * gd.d3m13 + 2.0*gd.qi2 * gd.qi3 * gd.d3m23)/2.0/gd.p0);

  
#if EVOLVE_LAMBDA
  p.x_c[6] = dt / gd.p0;
#endif
}

void Geodesic::registerRKRefiner(
  xfer::RefineAlgorithm& refiner,
  std::shared_ptr<hier::RefineOperator> &particle_refine_op)
{
  refiner.registerRefine(pc_idx,
                         pc_idx,                        
                         pc_idx,                        
                         particle_refine_op);

}

// set values needed by evolution equation for null geodesic
// since the corresponding particle has been living in certain
// patch, do not need global interpolation any more
void Geodesic::set_gd_values(
  const std::shared_ptr<hier::Patch> & patch, 
  double p_info[], GeodesicData *gd, BSSN *bssn, const real_t dx[], double shift[])
{
  gd->x = p_info[0], gd->y = p_info[1], gd->z = p_info[2];
  gd->q1 = p_info[3], gd->q2 = p_info[4], gd->q3 = p_info[5];
#if EVOLVE_LAMBDA
  gd->lambda = p_info[6];
#endif
  gd->x += (double)shift[0] * L[0];
  gd->y += (double)shift[1] * L[1];
  gd->z += (double)shift[2] * L[2];
  
  int i0 = floor((gd->x - domain_lower[0] ) / dx[0] - 0.5);
  int j0 = floor((gd->y - domain_lower[1] ) / dx[1] - 0.5);
  int k0 = floor((gd->z - domain_lower[2] ) / dx[2] - 0.5);

  
  
  real_t x0 = domain_lower[0] + (double)i0 * dx[0] + dx[0]/2.0;
  real_t y0 = domain_lower[1] + (double)j0 * dx[1] + dx[1]/2.0;
  real_t z0 = domain_lower[2] + (double)k0 * dx[2] + dx[2]/2.0;

  double xd = (gd->x - x0) / dx[0];
  double yd = (gd->y - y0) / dx[1];
  double zd = (gd->z - z0) / dx[2];

  GEODESIC_DEFINE_CRSPLINES_CHI;
  GEODESIC_DEFINE_CRSPLINES_DCHI(1);
  GEODESIC_DEFINE_CRSPLINES_DCHI(2);
  GEODESIC_DEFINE_CRSPLINES_DCHI(3);

  GEODESIC_DEFINE_CRSPLINES_BETA(1);
  GEODESIC_DEFINE_CRSPLINES_BETA(2);
  GEODESIC_DEFINE_CRSPLINES_BETA(3);

  GEODESIC_DEFINE_CRSPLINES_DBETA(1,1);
  GEODESIC_DEFINE_CRSPLINES_DBETA(1,2);
  GEODESIC_DEFINE_CRSPLINES_DBETA(1,3);
  GEODESIC_DEFINE_CRSPLINES_DBETA(2,1);
  GEODESIC_DEFINE_CRSPLINES_DBETA(2,2);
  GEODESIC_DEFINE_CRSPLINES_DBETA(2,3);
  GEODESIC_DEFINE_CRSPLINES_DBETA(3,1);
  GEODESIC_DEFINE_CRSPLINES_DBETA(3,2);
  GEODESIC_DEFINE_CRSPLINES_DBETA(3,3);

  GEODESIC_DEFINE_CRSPLINES_ALPHA;

  GEODESIC_DEFINE_CRSPLINES_DALPHA(1);
  GEODESIC_DEFINE_CRSPLINES_DALPHA(2);
  GEODESIC_DEFINE_CRSPLINES_DALPHA(3);

  COSMO_APPLY_TO_IJ_PERMS(GEODESIC_DEFINE_CRSPLINES_M);
  COSMO_APPLY_TO_IJK_PERMS(GEODESIC_DEFINE_CRSPLINES_DM);

        
  BSSNData bd = {0};

  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      for(int k = 0; k < 4; k++)
      {
#if USE_COSMOTRACE
       bssn->set_bd_values_for_ray_tracing(i0-1+i, j0-1+j, k0-1+k, &bd, dx);
#endif
        GEODESIC_CRSPLINES_SET_F_CHI;
        GEODESIC_CRSPLINES_SET_F_DCHI(1);
        GEODESIC_CRSPLINES_SET_F_DCHI(2);
        GEODESIC_CRSPLINES_SET_F_DCHI(3);

        GEODESIC_CRSPLINES_SET_F_BETA(1);
        GEODESIC_CRSPLINES_SET_F_BETA(2);
        GEODESIC_CRSPLINES_SET_F_BETA(3);

        GEODESIC_CRSPLINES_SET_F_DBETA(1,1);
        GEODESIC_CRSPLINES_SET_F_DBETA(1,2);
        GEODESIC_CRSPLINES_SET_F_DBETA(1,3);
        GEODESIC_CRSPLINES_SET_F_DBETA(2,1);
        GEODESIC_CRSPLINES_SET_F_DBETA(2,2);
        GEODESIC_CRSPLINES_SET_F_DBETA(2,3);
        GEODESIC_CRSPLINES_SET_F_DBETA(3,1);
        GEODESIC_CRSPLINES_SET_F_DBETA(3,2);
        GEODESIC_CRSPLINES_SET_F_DBETA(3,3);

        GEODESIC_CRSPLINES_SET_F_ALPHA;

        GEODESIC_CRSPLINES_SET_F_DALPHA(1);
        GEODESIC_CRSPLINES_SET_F_DALPHA(2);
        GEODESIC_CRSPLINES_SET_F_DALPHA(3);

        COSMO_APPLY_TO_IJ_PERMS(GEODESIC_CRSPLINES_SET_F_M);
        COSMO_APPLY_TO_IJK_PERMS(GEODESIC_CRSPLINES_SET_F_DM);
      }

  GEODESIC_CRSPLINES_CAL_COEF_CHI;
  GEODESIC_CRSPLINES_CAL_COEF_DCHI(1);
  GEODESIC_CRSPLINES_CAL_COEF_DCHI(2);
  GEODESIC_CRSPLINES_CAL_COEF_DCHI(3);

  GEODESIC_CRSPLINES_CAL_COEF_BETA(1);
  GEODESIC_CRSPLINES_CAL_COEF_BETA(2);
  GEODESIC_CRSPLINES_CAL_COEF_BETA(3);

  GEODESIC_CRSPLINES_CAL_COEF_DBETA(1,1);
  GEODESIC_CRSPLINES_CAL_COEF_DBETA(1,2);
  GEODESIC_CRSPLINES_CAL_COEF_DBETA(1,3);
  GEODESIC_CRSPLINES_CAL_COEF_DBETA(2,1);
  GEODESIC_CRSPLINES_CAL_COEF_DBETA(2,2);
  GEODESIC_CRSPLINES_CAL_COEF_DBETA(2,3);
  GEODESIC_CRSPLINES_CAL_COEF_DBETA(3,1);
  GEODESIC_CRSPLINES_CAL_COEF_DBETA(3,2);
  GEODESIC_CRSPLINES_CAL_COEF_DBETA(3,3);

  GEODESIC_CRSPLINES_CAL_COEF_ALPHA;

  GEODESIC_CRSPLINES_CAL_COEF_DALPHA(1);
  GEODESIC_CRSPLINES_CAL_COEF_DALPHA(2);
  GEODESIC_CRSPLINES_CAL_COEF_DALPHA(3);

  COSMO_APPLY_TO_IJ_PERMS(GEODESIC_CRSPLINES_CAL_COEF_M);
  COSMO_APPLY_TO_IJK_PERMS(GEODESIC_CRSPLINES_CAL_COEF_DM);


  
  
  GEODESIC_CRSPLINES_EVAL_CHI;
  GEODESIC_CRSPLINES_EVAL_DCHI(1);
  GEODESIC_CRSPLINES_EVAL_DCHI(2);
  GEODESIC_CRSPLINES_EVAL_DCHI(3);

  GEODESIC_CRSPLINES_EVAL_BETA(1);
  GEODESIC_CRSPLINES_EVAL_BETA(2);
  GEODESIC_CRSPLINES_EVAL_BETA(3);

  GEODESIC_CRSPLINES_EVAL_DBETA(1,1);
  GEODESIC_CRSPLINES_EVAL_DBETA(1,2);
  GEODESIC_CRSPLINES_EVAL_DBETA(1,3);
  GEODESIC_CRSPLINES_EVAL_DBETA(2,1);
  GEODESIC_CRSPLINES_EVAL_DBETA(2,2);
  GEODESIC_CRSPLINES_EVAL_DBETA(2,3);
  GEODESIC_CRSPLINES_EVAL_DBETA(3,1);
  GEODESIC_CRSPLINES_EVAL_DBETA(3,2);
  GEODESIC_CRSPLINES_EVAL_DBETA(3,3);

  GEODESIC_CRSPLINES_EVAL_ALPHA;

  GEODESIC_CRSPLINES_EVAL_DALPHA(1);
  GEODESIC_CRSPLINES_EVAL_DALPHA(2);
  GEODESIC_CRSPLINES_EVAL_DALPHA(3);

  COSMO_APPLY_TO_IJ_PERMS(GEODESIC_CRSPLINES_EVAL_M);
  COSMO_APPLY_TO_IJK_PERMS(GEODESIC_CRSPLINES_EVAL_DM);



  
  // calculating derivative to the 3-metric, not conformal one
  COSMO_APPLY_TO_IJK_PERMS(GEODESIC_CRSPLINES_CAL_DM);
  
  gd->m11 = gd->m11 / pw2(gd->chi);
  gd->m12 = gd->m12 / pw2(gd->chi);
  gd->m13 = gd->m13 / pw2(gd->chi);
  gd->m22 = gd->m22 / pw2(gd->chi);
  gd->m23 = gd->m23 / pw2(gd->chi);
  gd->m33 = gd->m33 / pw2(gd->chi);

  real_t det = gd->m11 * gd->m22 * gd->m33 + gd->m12 * gd->m23 * gd->m13
    + gd->m12 * gd->m23 * gd->m13 - gd->m13 * gd->m22 * gd->m13
    - gd->m12 * gd->m12 * gd->m33 - gd->m23 * gd->m23 * gd->m11;
  
  gd->mi11 = (gd->m22 * gd->m33 - pw2(gd->m23)) / det;
  gd->mi22 = (gd->m11 * gd->m33 - pw2(gd->m13)) / det;
  gd->mi33 = (gd->m11 * gd->m22 - pw2(gd->m12)) / det;
  gd->mi12 = (gd->m13*gd->m23 - gd->m12*(gd->m33)) / det;
  gd->mi13 = (gd->m12*gd->m23 - gd->m13*(gd->m22)) / det;
  gd->mi23 = (gd->m12*gd->m13 - gd->m23*(gd->m11)) / det;

  gd->qi1 = gd->q1 * gd->mi11 + gd->q2 * gd->mi12 + gd->q3 * gd->mi13;
  gd->qi2 = gd->q1 * gd->mi12 + gd->q2 * gd->mi22 + gd->q3 * gd->mi23;
  gd->qi3 = gd->q1 * gd->mi13 + gd->q2 * gd->mi23 + gd->q3 * gd->mi33;

  gd->p0 = sqrt(
    + gd->q1 * gd->q1 * gd->mi11 + 2.0*gd->q1 * gd->q2 * gd->mi12
    + gd->q2 * gd->q2 * gd->mi22 + 2.0*gd->q1 * gd->q3 * gd->mi13
    + gd->q3 * gd->q3 * gd->mi33 + 2.0*gd->q2 * gd->q3 * gd->mi23)/ gd->alpha;

  //  std::cout<<"!!! p0 "<<gd->p0<<" "<<gd->mi11<<" "<<gd->mi22<<" "<<gd->m11<<" "<<gd->chi<<"\n";
  
}

void Geodesic::K1FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);

  const hier::Box& box = patch->getBox();
  
  
  const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(  
    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  hier::Box ghost_box = pc_pdata->getGhostBox();

  pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iter(*pc_pdata, true);
  pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iterend(*pc_pdata, false);

  for(;iter != iterend; iter++)
  {

    ParticleContainer & id = *iter;

    for(std::list<RKParticle>::iterator it=id.p_list.begin();
        it != id.p_list.end(); it++)
    {
      for(int i = 0; i < PARTICLE_NUMBER_OF_STATES; i++)
      {
        (*it).x_a[i] = (*it).x_p[i] + (*it).x_c[i] / 2.0;
        (*it).x_f[i] = (*it).x_p[i] + (*it).x_c[i] / 6.0;
      }
      // if(patch->getBox().getBoxId().getOwnerRank() ==7)
      //   std::cout<<"haha "<<(*it).x_p[0]<<" "<<(*it).x_p[1]<<" "<<(*it).x_p[2]
      //            <<" "<<(*it).x_a[0]<<" "<<(*it).x_a[1]<<" "<<(*it).x_a[2]
      //            <<" "<<(*it).x_c[0]<<" "<<(*it).x_c[1]<<" "<<(*it).x_c[2]<<"\n"<<std::flush;

    }
    
  }
}

void Geodesic::K2FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);

  const hier::Box& box = patch->getBox();
  
  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
  const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(  
    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  //initialize dx for each patch
  const real_t * dx = &(patch_geom->getDx())[0];

  hier::Box ghost_box = pc_pdata->getGhostBox();


  pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iter(*pc_pdata, true);
  pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iterend(*pc_pdata, false);

  for(;iter != iterend; iter++)
  {

    ParticleContainer & id = *iter;

    for(std::list<RKParticle>::iterator it=id.p_list.begin();
        it != id.p_list.end(); it++)
    {
      for(int i = 0; i < PARTICLE_NUMBER_OF_STATES; i++)
      {
        (*it).x_a[i] = (*it).x_p[i] + (*it).x_c[i] / 2.0;
        (*it).x_f[i] += (*it).x_c[i] / 3.0;
      }
      
    }    
  }
}

void Geodesic::K3FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);

  const hier::Box& box = patch->getBox();
  
  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
  const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(  
    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  //initialize dx for each patch
  const real_t * dx = &(patch_geom->getDx())[0];

  hier::Box ghost_box = pc_pdata->getGhostBox();

  pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iter(*pc_pdata, true);
  pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iterend(*pc_pdata, false);

  for(;iter != iterend; iter++)
  {

    ParticleContainer & id = *iter;

    for(std::list<RKParticle>::iterator it=id.p_list.begin();
        it != id.p_list.end(); it++)
    {
      for(int i = 0; i < PARTICLE_NUMBER_OF_STATES; i++)
      {
        (*it).x_a[i] = (*it).x_p[i] + (*it).x_c[i];
        (*it).x_f[i] += (*it).x_c[i] / 3.0;
      }
      
    }

  }
}

void Geodesic::K4FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);

  const hier::Box& box = patch->getBox();
  
  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
  const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(  
    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  //initialize dx for each patch
  const real_t * dx = &(patch_geom->getDx())[0];

  hier::Box ghost_box = pc_pdata->getGhostBox();


  pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iter(*pc_pdata, true);
  pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iterend(*pc_pdata, false);

  for(;iter != iterend; iter++)
  {

    ParticleContainer &id = *iter;

    for(std::list<RKParticle>::iterator it=id.p_list.begin();
        it != id.p_list.end(); it++)
    {
      for(int i = 0; i < PARTICLE_NUMBER_OF_STATES; i++)
      {
        (*it).x_f[i] += (*it).x_c[i] / 6.0;
        (*it).x_a[i] = (*it).x_p[i] = (*it).x_f[i];
      }
      
    }

  }
}

void Geodesic::clearParticlesCoveredbyFinerLevel(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int ln)
{
  std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
  {
    std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));

    // loop over patches on level
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end(); ++ip)
    {
      std::shared_ptr<hier::Patch> patch(*ip);

      std::shared_ptr<pdat::CellData<double> > w_pdata(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_idx)));

      initPData(patch);
      
      std::shared_ptr<pdat::CellData<double> > w_cell_data_(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(w_pdata));

      pdat::CellData<double>& w_cell_data = *w_cell_data_;
       
      hier::Box::iterator iend(patch->getBox().end());

      for (hier::Box::iterator i(patch->getBox().begin()); i != iend; ++i)
      {
         const pdat::CellIndex cell_index(*i);
         //         tbox::pout<<"Cell "<<(*i)<<"has vol "<<w_cell_data(cell_index)<<"\n";
         ParticleContainer * p = pc_pdata->getItem(*i);
         
         if(p!= NULL && fabs(w_cell_data(cell_index)) < 1e-9) // it means it is covered by finer level
         {
           pc_pdata->removeItem(cell_index);
         }
      }

    }
  }

}


void Geodesic::regridPreProcessing(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  std::shared_ptr<hier::CoarsenOperator> &particle_coarsen_op)
{

  clearParticles(hierarchy, pc_s_idx);

  // coarsen to virtual particle containers
  xfer::CoarsenAlgorithm coarsener(dim);
    coarsener.registerCoarsen(pc_idx,                  
                              pc_s_idx,                   
                              particle_coarsen_op);

  for(int ln = 0; ln < hierarchy->getNumberOfLevels()-1; ln++)
  {
    std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));

    std::shared_ptr<xfer::CoarsenSchedule> coarsen_schedules =
      coarsener.createSchedule(level, hierarchy->getPatchLevel(ln+1));

    coarsen_schedules->coarsenData();
    
  }
}

void Geodesic::regridPostProcessing(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  for(int ln = 0; ln < hierarchy->getNumberOfLevels()-1; ln++)
  {
    // clear ghost particles and particles covered by finer grid
    clearParticlesCoveredbyFinerLevel(hierarchy, ln);
    // restore virtual particles those are NOT covered by finer grid anymore to
    // real particles
    std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));

        // loop over patches on level
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end(); ++ip)
    {
      std::shared_ptr<hier::Patch> patch(*ip);

      std::shared_ptr<pdat::CellData<double> > w_pdata(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_idx)));

      initPData(patch);
      
      std::shared_ptr<pdat::CellData<double> > w_cell_data_(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(w_pdata));

      pdat::CellData<double>& w_cell_data = *w_cell_data_;
       
      hier::Box::iterator iend(patch->getBox().end());

      for (hier::Box::iterator i(patch->getBox().begin()); i != iend; ++i)
      {
         const pdat::CellIndex cell_index(*i);
         //         tbox::pout<<"Cell "<<(*i)<<"has vol "<<w_cell_data(cell_index)<<"\n";
         ParticleContainer * src_p = pc_s_pdata->getItem(*i);


         if(src_p!= NULL && fabs(w_cell_data(cell_index)) > 1e-10) // it means it is NOT covered by finer level
         {

           ParticleContainer *dst_p = pc_pdata->getItem(*i);

           if(dst_p == NULL)
           {
             ParticleContainer *p = new ParticleContainer(src_p->idx);
             pc_pdata->appendItem((*i), *p);
             dst_p = pc_pdata->getItem(*i);
             
           }
           dst_p->p_list.insert(dst_p->p_list.end(), src_p->p_list.begin(), src_p->p_list.end());

         }
      }

    }
  }
  // clear virtual particles 
  clearParticles(hierarchy, pc_s_idx);

}

// particle redistribution, see the notes for more detail
void Geodesic::particleRedistribution(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  int ln, bool do_update_buffer, int step)
{
  //  if(step == 54) clearParticles(hierarchy, ln, pc_d_buffer_idx);

  cur_step = step;
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));
  // put all particles to correct cell after advancing it

  // first, step: clear scratch particles
  clearParticles(hierarchy, ln, pc_s_idx);
  
  for (hier::PatchLevel::iterator ip(level->begin());
       ip != level->end(); ++ip)
  {
    std::shared_ptr<hier::Patch> patch(*ip);
    initPData(patch);

    const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry()));

    const double *dx = &patch_geom->getDx()[0];

    const hier::Box& box = patch->getBox();
    hier::Box ghost_box = pc_pdata->getGhostBox();

      pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iter(*pc_pdata, true);
    pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iterend(*pc_pdata, false);
    
    // pdat::CellIterator icend(pdat::CellGeometry::end(ghost_box));
    // for (pdat::CellIterator ic(pdat::CellGeometry::begin(ghost_box));
    //      ic != icend; ++ic)
    for(; iter != iterend; iter++)
    {
      hier::Index old_idx(iter.getIndex());

      ParticleContainer & id = *iter;


      // go through all particles
      for(std::list<RKParticle>::iterator it=id.p_list.begin();
          it != id.p_list.end();)
      {
//         hier::Index temp_new_idx(floor((*it).x_a[0] / dx[0] ),
//                             floor((*it).x_a[1] / dx[1] ),
//                             floor((*it).x_a[2] / dx[2] ));

//         for(int i = 0; i < 3 ; i++)
//         {
//           if(temp_new_idx[i] > old_idx[i] + 1)
//             temp_new_idx[i] -= 
//               round((double)(temp_new_idx[i] - old_idx[i]) / round(L[i] / dx[i]))
//               * round(L[i] / dx[i]);
//           else if(temp_new_idx[i] < old_idx[i] - 1)
//             temp_new_idx[i] += 
//               round((double)(old_idx[i] - temp_new_idx[i]) / round(L[i] / dx[i]))
//               * round(L[i] / dx[i]);
//           if( abs(temp_new_idx[i] - old_idx[i]) > 1 )
//             TBOX_ERROR("The difference between old and new idx is larger than 1 "
// <<temp_new_idx<<" "<<old_idx<<"\n");
//         }

//         // still at the same position, nothing needs to be done
//         if(old_idx == temp_new_idx)
//         {
//           ++it;
//           continue;
//         }

        double shift[3] = {0};
        // shift new_idx for periodic boundary, VERY IMPORTANT!

        // calculating new cell index, should be very careful on periodic boundary
        // see following shift precedure 
        hier::Index new_idx(floor((*it).x_a[0] / dx[0] ),
                            floor((*it).x_a[1] / dx[1] ),
                            floor((*it).x_a[2] / dx[2] ));
        for(int i = 0; i < 3 ; i++)
        {
          if(new_idx[i] > old_idx[i] + 1)
            new_idx[i] -= 
              round((double)(new_idx[i] - old_idx[i]) / round(L[i] / dx[i]))
              * round(L[i] / dx[i]);
          else if(new_idx[i] < old_idx[i] - 1)
            new_idx[i] += 
              round((double)(old_idx[i] - new_idx[i]) / round(L[i] / dx[i]))
              * round(L[i] / dx[i]);
          if( abs(new_idx[i] - old_idx[i]) > 1 )
            TBOX_ERROR("The difference between old and new idx is larger than 1 "
                       <<new_idx<<" "<<old_idx<<"\n");
        }


        //std::cout<<mpi.getRank()<<" "<<new_idx<<" "<<old_idx<<" "<<(*it).x_a[2]<<"\n";
        // still at the same position, nothing needs to be done
        // if(old_idx == new_idx)
        // {
        //   ++it;
           
        //   continue;
        // }

        if(ghost_box.contains(new_idx))
        {
          ParticleContainer *fp = pc_s_pdata->getItem(new_idx);
          if(fp)
          {
            fp->addParticle(*it);
          }
          else
          {
            fp = new ParticleContainer(new_idx);
            fp->addParticle(*it);
            pc_s_pdata->appendItem(new_idx, *fp);
          }

        }
        else
        {
          // the new_idx is even out of the ghost region
          // theck whether they were ghost particles before advancing
          // if not, something went wrong
          if(box.contains(old_idx))
          {
            TBOX_ERROR("One particle was in real region but has gone out of the ghost region!!!\n");
          }
        }
                
        // std::cout<<"!!!!!\n";
        // std::cout<<(it==id.p_list.begin())<<" "<<(it_bak==id.p_list.begin())<<" "<<(id.p_list.size())<<"\n";
        // id.print();
        // erase and advance
        it = id.p_list.erase(it);
        
      }

      // if(id.p_list.empty())
      //   pc_pdata->removeItem(*ic);
    }
    pc_pdata->removeAllItems();
    insertPatchParticles(patch, pc_s_idx, pc_idx);
    pc_s_pdata->removeAllItems();
  }

  // updating upstream buffer
  if(ln < hierarchy->getNumberOfLevels() - 1)
  {
    std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));

    xfer::RefineAlgorithm refiner;

    std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
        hierarchy->getGridGeometry()));

    TBOX_ASSERT(grid_geometry_);

    geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;


    std::shared_ptr<hier::RefineOperator> refine_op
      = grid_geometry.
      lookupRefineOperator(pc, "PARTICLE_REFINE");

   
    refiner.registerRefine(pc_s_idx,
                           pc_idx,                        
                           pc_s_idx,                        
                           refine_op);

    std::shared_ptr<xfer::RefineSchedule> refine_schedule;

    std::shared_ptr<hier::PatchLevel> finer_level(
      hierarchy->getPatchLevel(ln + 1));

    // fill only real cells of finer level
    refine_schedule =
      refiner.createSchedule(
        std::shared_ptr<xfer::PatchLevelFillPattern>(
          new xfer::PatchLevelInteriorFillPattern()),
        finer_level,
        NULL,
        ln,
        hierarchy,
        NULL);


    refine_schedule->fillData(0.0);

    /*********************************************/
  

    // xfer::RefineAlgorithm refiner2;
  
    // refiner2.registerRefine(scratch_id,
    //                        scratch_id,                        
    //                        scratch_id,                        
    //                        refine_op);

    // std::shared_ptr<xfer::RefineSchedule> refine_schedule2;
  
    // refine_schedule2 =
    //   refiner.createSchedule(
    //     level,
    //     NULL);

    // refine_schedule2->fillData(0.0);

    /*****************************************/    
    insertLevelParticles(hierarchy, ln + 1, pc_s_idx, pc_idx);

    // clear scratch particles
    clearParticles(hierarchy, ln + 1, pc_s_idx);
    
    clearParticlesCoveredbyFinerLevel(hierarchy, ln);

    
    //coarsen the data from downstream buffer from the upper level
    xfer::CoarsenAlgorithm coarsener(dim);

    std::shared_ptr<hier::CoarsenOperator> coarsen_op = grid_geometry.
      lookupCoarsenOperator(pc, "PARTICLE_COARSEN");

    coarsener.registerCoarsen(pc_s_idx,                  
                              pc_d_buffer_idx,                   
                              coarsen_op,
                              hier::IntVector(dim, 1));

    std::shared_ptr<xfer::CoarsenSchedule> coarsen_schedules =
      coarsener.createSchedule(level, hierarchy->getPatchLevel(ln+1));

    coarsen_schedules->coarsenData();

    insertLevelParticles(hierarchy, ln, pc_s_idx, pc_idx);

    //    coarsen_schedules->printClassData(tbox::plog);
    
    // clear scatch particles
    clearParticles(hierarchy, ln, pc_s_idx);
  }

  // updating downstream buffer
  if(ln > 0 && do_update_buffer)
  {

        int p_in_buffer_num = 0;
    std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));
    //    if(step == 54) std::cout<<"Flag2.0000 "<<pc_d_buffer_idx<<"\n"<<std::flush;      

    clearParticles(hierarchy, ln, pc_d_buffer_idx);


    //  if(step == 54) std::cout<<"Flag2.0!"<<std::flush;      
    // the CoarseFineBoundary type deals with coarse-fine boundary boxes on certain level
    hier::CoarseFineBoundary cfb(*hierarchy, ln, hier::IntVector(dim, ghost_width));       
    // cfb.printClassData(tbox::pout);
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end(); ++ip)
    {
      std::shared_ptr<hier::Patch> patch(*ip);


      
      std::shared_ptr<hier::PatchGeometry> geom (patch->getPatchGeometry());


      initPData(patch);

      // enumerating boudary boxes with different co-dimension as coarse-fine boundaries
      
      const std::vector<hier::BoundaryBox> & codim1_boxes =
        cfb.getBoundaries(patch->getGlobalId(), 1, patch->getBox().getBlockId());

      //  if(step == 54) std::cout<<"Flag2.1!"<<std::flush;  
      const int n_codim1_boxes = static_cast<int>(codim1_boxes.size());

      for(int l = 0 ; l < n_codim1_boxes; l++)
      {
        hier::Box boundary_fill_box =
          geom->getBoundaryFillBox(
            codim1_boxes[l], patch->getBox(), hier::IntVector(dim, ghost_width));

        pdat::CellIterator icend(pdat::CellGeometry::end(boundary_fill_box));
        for (pdat::CellIterator ic(pdat::CellGeometry::begin(boundary_fill_box));
             ic != icend; ++ic)
        {
          ParticleContainer *pc = pc_pdata->getItem(*ic);
          if(pc != NULL){
            pc_d_buffer_pdata->replaceAddItem(*ic, *pc);

          }
        }
      }
      //  if(step == 54) std::cout<<"Flag2.2!"<<std::flush;  
      const std::vector<hier::BoundaryBox> & codim2_boxes =
        cfb.getBoundaries(patch->getGlobalId(), 2, patch->getBox().getBlockId());


      const int n_codim2_boxes = static_cast<int>(codim2_boxes.size());

      for(int l = 0 ; l < n_codim2_boxes; l++)
      {
        hier::Box boundary_fill_box =
          geom->getBoundaryFillBox(
            codim2_boxes[l], patch->getBox(), hier::IntVector(dim, ghost_width));


        pdat::CellIterator icend(pdat::CellGeometry::end(boundary_fill_box));
        for (pdat::CellIterator ic(pdat::CellGeometry::begin(boundary_fill_box));
             ic != icend; ++ic)
        {
          ParticleContainer *pc = pc_pdata->getItem(*ic);
          if(pc != NULL)
          {
            //            std::cout<<"Inserting "<<step<<" "<<(*ic)<<"\n"<<std::flush;

            pc_d_buffer_pdata->replaceAddItem(*ic, *pc);
          }
        }
        
        
      }

      //  if(step == 54) std::cout<<"Flag2.3!"<<std::flush;  
      const std::vector<hier::BoundaryBox> & codim3_boxes =
        cfb.getBoundaries(patch->getGlobalId(), 3, patch->getBox().getBlockId());


      const int n_codim3_boxes = static_cast<int>(codim3_boxes.size());

      for(int l = 0 ; l < n_codim3_boxes; l++)
      {
        hier::Box boundary_fill_box =
          geom->getBoundaryFillBox(
            codim3_boxes[l], patch->getBox(), hier::IntVector(dim, ghost_width));

        pdat::CellIterator icend(pdat::CellGeometry::end(boundary_fill_box));
        for (pdat::CellIterator ic(pdat::CellGeometry::begin(boundary_fill_box));
             ic != icend; ++ic)
        {
          ParticleContainer *pc = pc_pdata->getItem(*ic);

          if(pc != NULL){
            pc_d_buffer_pdata->replaceAddItem(*ic, *pc);
          }
        }

                
      }
      p_in_buffer_num += pc_d_buffer_pdata->getNumberOfItems();

    }

    if(mpi.getSize() > 1)
      mpi.AllReduce(&p_in_buffer_num, 1, MPI_SUM);

    //    std::cout<<"Number of particles in buffer is "<<p_in_buffer_num<<"\n"<<std::flush;

  }

}

// print particles info with certain id 
void Geodesic::printAll(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int idx)
{
  tbox::pout<<"Going to print all particles info in hierarchy\n";

  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());

  for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
  {
    int cnt = 0;
    
    tbox::pout<<"Starting to print particles on level "<<ln<<"\n";
    std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));

    
    // loop over patches on level
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end(); ++ip)
    {
      std::shared_ptr<hier::Patch> patch(*ip);
    
      // access sample data from patch
      std::shared_ptr<pdat::IndexData<ParticleContainer,
                                        pdat::CellGeometry> > sample(
                                          SAMRAI_SHARED_PTR_CAST<pdat::IndexData<ParticleContainer, pdat::CellGeometry>,
                                          hier::PatchData>(
                                            patch->getPatchData(idx)));
      TBOX_ASSERT(sample);


      // iterate over the index data stored on the patch
      // and dump the particles stored on it.

      pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iter(*sample, true);
      pdat::IndexData<ParticleContainer, pdat::CellGeometry>::iterator iterend(*sample, false);

      // pdat::CellIterator icend(pdat::CellGeometry::end(sample->getGhostBox()));
      // for (pdat::CellIterator ic(pdat::CellGeometry::begin(sample->getGhostBox()));
      //      ic != icend; ++ic) {
      for(; iter != iterend; iter++){
        ParticleContainer & p = *iter;
        {
          std::cout<<"\nThere are "<<p.p_list.size()<<" particles at cell "
                   <<iter.getIndex()<<". Node "<<mpi.getRank()<<", box "<<patch->getBox()
                   <<". Their location and velocity are \n";
          cnt += p.p_list.size();
          p.print();
        }
      }
    }
    if (mpi.getSize() > 1)
    {
      mpi.AllReduce(&cnt, 1, MPI_SUM);
    }

    tbox::pout<<"Finish outputing "<<cnt<<" particles on level "<<ln<<"\n";

  }

}


}
