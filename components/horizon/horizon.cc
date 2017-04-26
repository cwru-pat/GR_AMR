#include "horizon.h"
#include "../../cosmo_includes.h"
#include "../../utils/math.h"
#include "../bssn/bssn_data.h"

using namespace SAMRAI;

namespace cosmo
{
  
Horizon::Horizon(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  const tbox::Dimension& dim_in,
  boost::shared_ptr<tbox::Database> database_in,
  std::ostream* l_stream_in):
  lstream(l_stream_in),
  cosmo_horizon_db(database_in),
  dim(dim_in),
  is_periodic(cosmo_horizon_db->getBoolWithDefault("is_periodic", false))
{
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

  boost::shared_ptr<hier::VariableContext> context_scratch(
    variable_db->getContext("SCRATCH"));
  boost::shared_ptr<hier::VariableContext> context_active(
    variable_db->getContext("ACTIVE"));
  boost::shared_ptr<hier::VariableContext> context_previous(
    variable_db->getContext("PREVIOUS"));
  boost::shared_ptr<hier::VariableContext> context_k1(
    variable_db->getContext("RK_K1"));
  boost::shared_ptr<hier::VariableContext> context_k2(
    variable_db->getContext("RK_K2"));
  boost::shared_ptr<hier::VariableContext> context_k3(
    variable_db->getContext("RK_K3"));
  boost::shared_ptr<hier::VariableContext> context_k4(
    variable_db->getContext("RK_K4"));

  
  VAR_INIT(F);
  VAR_INIT(s1);
  VAR_INIT(s2);
  VAR_INIT(s3);

  
  REG_TO_CONTEXT(F, context_scratch, s, STENCIL_ORDER);
  REG_TO_CONTEXT(F, context_previous, p, STENCIL_ORDER);
  REG_TO_CONTEXT(F, context_active, a, STENCIL_ORDER);
  REG_TO_CONTEXT(F, context_k1, k1, STENCIL_ORDER);
  REG_TO_CONTEXT(F, context_k2, k2, STENCIL_ORDER);
  REG_TO_CONTEXT(F, context_k3, k3, STENCIL_ORDER);
  REG_TO_CONTEXT(F, context_k4, k4, STENCIL_ORDER);

  REG_TO_CONTEXT(s1, context_active, a, STENCIL_ORDER);
  REG_TO_CONTEXT(s2, context_active, a, STENCIL_ORDER);
  REG_TO_CONTEXT(s3, context_active, a, STENCIL_ORDER);

  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
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
  origin = cosmo_horizon_db->getDoubleVector("origin");
  radius = cosmo_horizon_db->getDouble("radius");

}
Horizon::~Horizon()
{
  
}

void Horizon::clear(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln)
{
  math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);
  RK4_ARRAY_ZERO(F, hcellmath);
  EXTRA_ARRAY_ZERO(s2, hcellmath);
  EXTRA_ARRAY_ZERO(s2, hcellmath);
  EXTRA_ARRAY_ZERO(s3, hcellmath);
}

  
void Horizon::initSurface(const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{

  int ln_num = hierarchy->getNumberOfLevels();
  if(is_periodic)
  {
    for (int ln = 0; ln < ln_num; ln++)
    {
      boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    
      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
        const boost::shared_ptr<hier::Patch>& patch = *p;
        boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
          BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch->getPatchGeometry()));


        initPData(patch);
        initMDA(patch);
        const hier::Box& box = F_a_pdata->getGhostBox();
  
        const int * lower = &box.lower()[0];
        const int * upper = &box.upper()[0];

        const real_t * dx = &(patch_geometry->getDx())[0];

     
      
        for(int k = lower[2]; k <= upper[2]; k++)
        {
          for(int j = lower[1]; j <= upper[1]; j++)
          {
            for(int i = lower[0]; i <= upper[0]; i++)
            {
              real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0;
              real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0;
              real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0;
              real_t min_r = pw2(x - origin[0]) + pw2(y - origin[1]) + pw2(z - origin[2]);
              min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_lower[0]) + pw2(y - domain_lower[1]) + pw2(z - domain_lower[2]));
              min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_upper[0]) + pw2(y - domain_lower[1]) + pw2(z - domain_lower[2]));
              min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_lower[0]) + pw2(y - domain_upper[1]) + pw2(z - domain_lower[2]));
              min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_lower[0]) + pw2(y - domain_lower[1]) + pw2(z - domain_upper[2]));
              min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_lower[0]) + pw2(y - domain_upper[1]) + pw2(z - domain_upper[2]));
              min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_upper[0]) + pw2(y - domain_lower[1]) + pw2(z - domain_upper[2]));
              min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_upper[0]) + pw2(y - domain_upper[1]) + pw2(z - domain_lower[2]));
              min_r = tbox::MathUtilities<double>::Min(min_r, pw2(x - domain_upper[0]) + pw2(y - domain_upper[1]) + pw2(z - domain_upper[2]));

              min_r = sqrt(min_r);
            
              F_a(i, j, k) = F_p(i, j, k) =
                min_r - radius; 
            }
          }
        }     
      }

    }  // loop over levels

    return;
  }
  for (int ln = 0; ln < ln_num; ln++)
  {
    boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    
    for (hier::PatchLevel::iterator p(level->begin());
         p != level->end(); ++p) {
      const boost::shared_ptr<hier::Patch>& patch = *p;
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));


      initPData(patch);
      initMDA(patch);
      const hier::Box& box = F_a_pdata->getGhostBox();
  
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];

      const real_t * dx = &(patch_geometry->getDx())[0];

     
      
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0;
            real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0;
            real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0;
            F_a(i, j, k) = F_p(i, j, k) =
              sqrt(pw2(x - origin[0]) + pw2(y - origin[1]) + pw2(z - origin[2])) - radius; 
          }
        }
      }     
    }

  }  // loop over levels

}

/**
 * @brief initilizing all pointers for every components's patchdata
 * 
 */
void Horizon::initPData(
  const boost::shared_ptr<hier::Patch> & patch)
{
  PDATA_ALL_INIT(F);
  PDATA_INIT(s1, a);
  PDATA_INIT(s2, a);
  PDATA_INIT(s3, a);
}


/**
 * @brief initilizing all pointers for every components's array access
 * 
 */
void Horizon::initMDA(
  const boost::shared_ptr<hier::Patch> & patch)
{

  MDA_ACCESS_ALL_INIT(F);
  MDA_ACCESS_INIT(s1, a);
  MDA_ACCESS_INIT(s2, a);
  MDA_ACCESS_INIT(s3, a);
}

void Horizon::setLevelTime(
  const boost::shared_ptr<hier::PatchLevel> & level,
  double from_t, double to_t)
{
  SET_LEVEL_TIME(F, from_t, to_t);
}

  
HorizonData Horizon::getHorizonData(
  idx_t i, idx_t j, idx_t k, BSSNData *bd, const real_t dx[])
{
  HorizonData hd = {0};

  hd.F = F_a(i, j, k);
  hd.s1 = s1_a(i, j, k);
  hd.s2 = s2_a(i, j, k);
  hd.s3 = s3_a(i, j, k);

  hd.d1F = derivative(i, j, k, 1, F_a, dx);
  hd.d2F = derivative(i, j, k, 2, F_a, dx);
  hd.d3F = derivative(i, j, k, 3, F_a, dx);

  hd.d1d1F = double_derivative(i, j, k, 1, 1, F_a, dx);
  hd.d2d2F = double_derivative(i, j, k, 2, 2, F_a, dx);
  hd.d3d3F = double_derivative(i, j, k, 3, 3, F_a, dx);
  hd.d1d2F = double_derivative(i, j, k, 1, 2, F_a, dx);
  hd.d1d3F = double_derivative(i, j, k, 1, 3, F_a, dx);
  hd.d2d3F = double_derivative(i, j, k, 2, 3, F_a, dx);
  
  return hd;
}

void Horizon::addNormVector(
  BSSN * bssn, const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  for(idx_t ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
  {
    const boost::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));
    addNormVector(bssn, level);
  }
  
}

void Horizon::addNormVector(
  BSSN * bssn, const boost::shared_ptr<hier::PatchLevel> & level)
{
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    addNormVector(bssn, patch);
  }  
}

  
void Horizon::addNormVector(
  BSSN * bssn, const boost::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);
  bssn->initPData(patch);
  bssn->initMDA(patch);

  
  const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom( 
    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  const real_t * dx = &(patch_geom->getDx())[0];

  const hier::Box& box = patch->getBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        BSSNData bd = {0};
        // TODO: remove redundant computations here?
        bssn->set_bd_values(i, j, k, &bd, dx);
        HorizonData hd = getHorizonData(i, j, k, &bd, dx);

        s1_a(i, j, k) = bd.chi *
          (bd.gammai11 * hd.d1F + bd.gammai12 * hd.d2F + bd.gammai13 * hd.d3F) /
          sqrt( 
               (bd.gammai11 * hd.d1F * hd.d1F + bd.gammai22 * hd.d2F * hd.d2F + bd.gammai33 * hd.d3F *hd.d3F
                + 2.0 * (bd.gammai12 * hd.d1F * hd.d2F + bd.gammai13 * hd.d1F * hd.d3F + bd.gammai23 * hd.d2F * hd.d3F)));

        s2_a(i, j, k) = bd.chi *
          (bd.gammai21 * hd.d1F + bd.gammai22 * hd.d2F + bd.gammai23 * hd.d3F) /
          sqrt(
               (bd.gammai11 * hd.d1F * hd.d1F + bd.gammai22 * hd.d2F * hd.d2F + bd.gammai33 * hd.d3F *hd.d3F
                + 2.0 * (bd.gammai12 * hd.d1F * hd.d2F + bd.gammai13 * hd.d1F * hd.d3F + bd.gammai23 * hd.d2F * hd.d3F)));
        s3_a(i, j, k) = bd.chi *
          (bd.gammai31 * hd.d1F + bd.gammai32 * hd.d2F + bd.gammai33 * hd.d3F) /
          sqrt(
               (bd.gammai11 * hd.d1F * hd.d1F + bd.gammai22 * hd.d2F * hd.d2F + bd.gammai33 * hd.d3F *hd.d3F
                + 2.0 * (bd.gammai12 * hd.d1F * hd.d2F + bd.gammai13 * hd.d1F * hd.d3F + bd.gammai23 * hd.d2F * hd.d3F)));

      }
    }
  }
  return;      
}

real_t Horizon::ev_F(BSSNData *bd, HorizonData *hd, const real_t dx[])
{
  real_t diFdiF = pw2(bd->chi) * (hd->d1F * hd->d1F * bd->gammai11 + hd->d2F * hd->d2F * bd->gammai22 + hd->d3F * hd->d3F * bd->gammai33
         + 2.0 * (hd->d1F * hd->d2F * bd->gammai12 + hd->d1F * hd->d3F * bd->gammai13 + hd->d2F * hd->d3F * bd->gammai23 ));

  
  if(bd->i == 32 && bd->j == 32 && bd->k ==32)
  {

    
    std::cout<<((pw2(bd->chi) * COSMO_SUMMATION_2(HORIZON_CALCULATE_DS0)) / sqrt(diFdiF)
           - pw2(bd->chi) * pw2(bd->chi) * COSMO_SUMMATION_2(HORIZON_CALCULATE_DS) / pow(diFdiF, 1.5)
    - bd->K
    +  1.0 / pw2(bd->chi) * (
      (bd->A11 + bd->K * bd->gamma11 / 3.0) * hd->s1 * hd->s1
      + (bd->A22 + bd->K * bd->gamma22 / 3.0) * hd->s2 * hd->s2
      + (bd->A33 + bd->K * bd->gamma33 / 3.0) * hd->s3 * hd->s3
      + 2.0 * (bd->A12 + bd->K * bd->gamma12 / 3.0) * hd->s1 * hd->s2
      + 2.0 * (bd->A13 + bd->K * bd->gamma13 / 3.0) * hd->s1 * hd->s3
      + 2.0 * (bd->A23 + bd->K * bd->gamma23 / 3.0) * hd->s2 * hd->s3
    ))<<"\n";
  }
                
  
  return ((pw2(bd->chi) * COSMO_SUMMATION_2(HORIZON_CALCULATE_DS0)) / sqrt(diFdiF)
           - pw2(bd->chi) * pw2(bd->chi) * COSMO_SUMMATION_2(HORIZON_CALCULATE_DS) / pow(diFdiF, 1.5)
    - bd->K
    +  1.0 / pw2(bd->chi) * (
      (bd->A11 + bd->K * bd->gamma11 / 3.0) * hd->s1 * hd->s1
      + (bd->A22 + bd->K * bd->gamma22 / 3.0) * hd->s2 * hd->s2
      + (bd->A33 + bd->K * bd->gamma33 / 3.0) * hd->s3 * hd->s3
      + 2.0 * (bd->A12 + bd->K * bd->gamma12 / 3.0) * hd->s1 * hd->s2
      + 2.0 * (bd->A13 + bd->K * bd->gamma13 / 3.0) * hd->s1 * hd->s3
      + 2.0 * (bd->A23 + bd->K * bd->gamma23 / 3.0) * hd->s2 * hd->s3
    )) * bd->chi *
    sqrt(hd->d1F * hd->d1F * bd->gammai11 + hd->d2F * hd->d2F * bd->gammai22 + hd->d3F * hd->d3F * bd->gammai33
         + 2.0 * (hd->d1F * hd->d2F * bd->gammai12 + hd->d1F * hd->d3F * bd->gammai13 + hd->d2F * hd->d3F * bd->gammai23 ))  ;
}    

void Horizon::RKEvolvePt(
  idx_t i, idx_t j, idx_t k, BSSNData &bd, const real_t dx[], real_t dt)
{
  HorizonData hd = getHorizonData(i, j, k, &bd, dx);

  F_s(i,j,k) = ev_F(&bd, &hd, dx) * dt;
}

bool Horizon::onTheSurface(idx_t i, idx_t j, idx_t k)
{
  if(SIGN(F_a(i-1, j, k) * F_a(i+1, j, k)) < 0 
     || SIGN(F_a(i, j-1, k) * F_a(i, j+1, k)) < 0
     || SIGN(F_a(i, j, k-1) * F_a(i, j, k+1)) < 0)
    return true;
  
  return false;
}

void Horizon::prepareForK1(
  const boost::shared_ptr<hier::PatchLevel> & level,
  real_t to_t)
{
  if(level == NULL) return;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    
    initPData(patch);
    initMDA(patch);

    const hier::Box& box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];


    // right branch of the tree
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(F_a_idx)->getTime()) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_R_K1(F);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(F_p_idx)->getTime())
       - (patch->getPatchData(F_a_idx)->getTime()
          - patch->getPatchData(F_p_idx)->getTime())) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_L_K1(F);
          }
        }
      }
    }
    else
      TBOX_ERROR("Current level locates neigher L or R branch of its father, check your code!");
    
  }
}

/**
 * @brief calculate K2 value on coarser level
 * 
 */
void Horizon::prepareForK2(
  const boost::shared_ptr<hier::PatchLevel> & level,
  real_t to_t)
{
  if(level == NULL) return;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    initPData(patch);
    initMDA(patch);

    const hier::Box& box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];


    // right branch of the tree
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(F_a_idx)->getTime()) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_R_K2(F);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(F_p_idx)->getTime())
       - (patch->getPatchData(F_a_idx)->getTime()
          - patch->getPatchData(F_p_idx)->getTime())) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_L_K2(F);
          }
        }
      }
    }
  }
}

/**
 * @brief calculate K3 value on coarser level
 * 
 */  
void Horizon::prepareForK3(
  const boost::shared_ptr<hier::PatchLevel> & level,
  real_t to_t)
{
  if(level == NULL) return;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    initPData(patch);
    initMDA(patch);

    const hier::Box& box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    // right branch of the tree
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(F_a_idx)->getTime()) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_R_K3(F);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(F_p_idx)->getTime())
       - (patch->getPatchData(F_a_idx)->getTime()
          - patch->getPatchData(F_p_idx)->getTime())) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_L_K3(F);
          }
        }
      }
    }
  }
}

/**
 * @brief calculate K4 value on coarser level
 * 
 */
void Horizon::prepareForK4(
  const boost::shared_ptr<hier::PatchLevel> & level,
  double to_t)
{
  if(level == NULL) return;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    initPData(patch);
    initMDA(patch);

    const hier::Box& box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    
    


    // right branch of the tree
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(F_a_idx)->getTime()) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_R_K4(F);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(F_p_idx)->getTime())
       - (patch->getPatchData(F_a_idx)->getTime()
          - patch->getPatchData(F_p_idx)->getTime())) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            RK4_INIT_L_K4(F);
          }
        }
      }
    }
  }
}  

void Horizon::alloc(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln)
{
  boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
  RK4_ARRAY_ALLOC(F);
  EXTRA_ARRAY_ALLOC(s1);
  EXTRA_ARRAY_ALLOC(s2);
  EXTRA_ARRAY_ALLOC(s3);
}


void Horizon::registerRKRefinerActive(
  xfer::RefineAlgorithm& refiner,
  boost::shared_ptr<hier::RefineOperator> &space_refine_op)
{
  REGISTER_SPACE_REFINE_A(F, refiner, space_refine_op);
}

  
/**
 * @brief register "scratch" component of the fields to refiner
 * 
 * @param refiner
 * @param refiner operator
 */
void Horizon::registerRKRefiner(
  xfer::RefineAlgorithm& refiner,
  boost::shared_ptr<hier::RefineOperator> &space_refine_op)
{
  REGISTER_SPACE_REFINE_S(F, refiner, space_refine_op);
}

/**
 * @brief register "active" component of the fields to coarsener 
 * 
 * @param coarsener
 * @param coarse operator
 */
void Horizon::registerCoarsenActive(
  xfer::CoarsenAlgorithm& coarsener,
  boost::shared_ptr<hier::CoarsenOperator>& coarsen_op)
{
  REGISTER_COARSEN_A(F, coarsener, coarsen_op);
}


void Horizon::copyAToP(
  math::HierarchyCellDataOpsReal<real_t> & hcellmath)
{
  COPY_A_TO_P(F);
}

/**
 * @brief  finalize k1 step for RK on both interior and boundary
 * 
 */
void Horizon::K1FinalizePatch(
  const boost::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);
  
  const hier::Box& box = F_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        RK4_FINALIZE_FIELD_1(F);
      }
    }
  }
}

void Horizon::K2FinalizePatch(
  const boost::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  const hier::Box& box = F_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        RK4_FINALIZE_FIELD_2(F);
      }
    }
  }
}


void Horizon::K3FinalizePatch(
  const boost::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  const hier::Box& box = F_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        RK4_FINALIZE_FIELD_3(F);
      }
    }
  }
}


void Horizon::K4FinalizePatch(
  const boost::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  //Only evolv inner grids for K4
  const hier::Box& box = F_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        RK4_FINALIZE_FIELD_4(F);
      }
    }
  }
}


void Horizon::RKEvolveHorizon(
  const boost::shared_ptr<hier::Patch> & patch, BSSN * bssn, real_t dt)
{
  bssn->initPData(patch);
  bssn->initMDA(patch);
  initPData(patch);
  initMDA(patch);

  
  const hier::Box& box = patch->getBox();
  
  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  BSSNData bd = {0};
  
  const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(  
    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  //initialize dx for each patch
  const real_t * dx = &(patch_geom->getDx())[0];


  
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        if( i == 32 && j == 32 && k == 32)
          std::cout<<F_a(i, j, k)<<" ";
        bssn->set_bd_values(i, j, k, &bd, dx);
        RKEvolvePt(i, j, k, bd, dx, dt);
        
      }
    }
  }
}

real_t Horizon::maxSurfaceMove(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln, idx_t w_idx)
{
  real_t max_change = -1;
  boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

  
  for (hier::PatchLevel::iterator p(level->begin());
       p != level->end(); ++p)
  {
    const boost::shared_ptr<hier::Patch>& patch = *p;
    boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
      BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry()));
  
    initPData(patch);

    initMDA(patch);

    boost::shared_ptr<pdat::CellData<real_t>> w_pdata (    
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  
        patch->getPatchData(w_idx)));
    arr_t w = pdat::ArrayDataAccess::access<DIM, double>(  
      w_pdata->getArrayData());
    
    const hier::Box& box = patch->getBox();
  
    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    const real_t * dx = &(patch_geometry->getDx())[0];

    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          if(w(i, j, k) > 0 && onTheSurface(i, j, k)
             && fabs(F_a(i, j, k) - F_p(i, j, k)) > max_change)
            max_change = fabs(F_a(i, j, k) - F_p(i, j, k));
        }
      }
    }


  }
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  mpi.Barrier();
  if (mpi.getSize() > 1) {
    mpi.AllReduce(&max_change, 1, MPI_MAX);
  }

  return max_change;
  
}

real_t Horizon::getFonBD(
  const boost::shared_ptr<hier::Patch> & patch,
  real_t bx, real_t by, real_t bz,
  hier::Box& g_box)
{
  initPData(patch);
  initMDA(patch);

  const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom( 
    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  const real_t * dx = &(patch_geom->getDx())[0];

  
  const int * lower = &g_box.lower()[0];
  const int * upper = &g_box.upper()[0];

  real_t min_s = INF;
  idx_t min_i = lower[0], min_j = lower[1], min_k = lower[2];
  real_t min_r = -1.0;
  
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0;
        real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0;
        real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0;

        real_t a = (x - origin[0]) * (by - origin[1]) - (y - origin[1]) * (bx - origin[0]);
        real_t b = (x - origin[0]) * (bz - origin[2]) - (z - origin[2]) * (bx - origin[0]);
        real_t c = (y - origin[1]) * (bz - origin[2]) - (z - origin[2]) * (by - origin[1]);
        
        real_t  slope = tbox::MathUtilities<double>::Max(
          tbox::MathUtilities<double>::Max(
            tbox::MathUtilities<double>::Abs(a), tbox::MathUtilities<double>::Abs(b)),
          tbox::MathUtilities<double>::Abs(c));
        if(slope < min_s)
        {
          min_s = slope;
          min_i = i, min_j = j, min_k = k;
          min_r = sqrt(pw2(x - origin[0]) + pw2(y - origin[1]) + pw2(z - origin[2]));
        }
      }
    }
  }

  if(min_r < 0)
    TBOX_ERROR("Cannot find grid with right slope!");

  return min_r - F_a(min_i, min_j, min_k);
  
}
  
void Horizon::updateBD(
  const boost::shared_ptr<hier::Patch> & patch,
  real_t dt)
{
  boost::shared_ptr<hier::PatchGeometry> geom (patch->getPatchGeometry());

  idx_t codim = 1;

  // getting all codimension 1 boxes
  const std::vector<hier::BoundaryBox> & codim1_boxes =
    geom->getCodimensionBoundaries(codim);

  const idx_t n_codim1_boxes = static_cast<idx_t>(codim1_boxes.size());

  // if it has no codimention 1 boundary, it has no other type of boundaries
  if(n_codim1_boxes == 0) return;

  initPData(patch);
  initMDA(patch);

  hier::Box& gg_box = const_cast<hier::Box&> (F_a_pdata->getGhostBox());

  hier::Box g_box = gg_box;
  
  const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom( 
    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

    
  const hier::Box& patch_box = patch->getBox();

  for(int i = 0 ; i < n_codim1_boxes; i++)
  {
    hier::Box boundary_fill_box =
      geom->getBoundaryFillBox(
        codim1_boxes[i], patch_box, F_a_pdata->getGhostCellWidth());

    if(boundary_fill_box.empty()) continue;


    idx_t l_idx = codim1_boxes[i].getLocationIndex();

    if(l_idx%2)
      g_box.growUpper(
        (hier::Box::dir_t)l_idx/2, -2*STENCIL_ORDER_WIDTH);
    else
      g_box.growLower(
        (hier::Box::dir_t)l_idx/2, -2*STENCIL_ORDER_WIDTH);
  }

  
  for(int l = 0 ; l < n_codim1_boxes; l++)
  {
    hier::Box boundary_fill_box =
      geom->getBoundaryFillBox(
        codim1_boxes[l], patch_box, F_a_pdata->getGhostCellWidth());

    if(boundary_fill_box.empty()) continue;


    idx_t l_idx = codim1_boxes[l].getLocationIndex();
    
    BSSNData bd = {0};

    boundary_fill_box.shift(
      (hier::Box::dir_t)l_idx/2,
      (l_idx%2)?(-STENCIL_ORDER_WIDTH):STENCIL_ORDER_WIDTH);


    boundary_fill_box *= patch_box;

    //initialize dx for each patch
    const real_t * dx = &(patch_geom->getDx())[0];

    const idx_t * lower = &boundary_fill_box.lower()[0];
    const idx_t * upper = &boundary_fill_box.upper()[0];
    if(l_idx == 0)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = upper[0]; i >= lower[0]; i--)
          {
            F_a(i, j, k) = 3.0 * (F_a(i+1, j, k) - F_a(i+2, j, k)) + F_a(i+3, j, k);

          }
        }
      }
    }
    if(l_idx == 1)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            F_a(i, j, k) = 3.0 * (F_a(i-1, j, k) - F_a(i-2, j, k)) + F_a(i-3, j, k);

          }
        }
      }
    }
    if(l_idx == 2)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = upper[1]; j >= lower[1]; j--)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            F_a(i, j, k) = 3.0 * (F_a(i, j+1, k) - F_a(i, j+2, k)) + F_a(i, j+3, k);

          }
        }
      }
    }
    if(l_idx == 3)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            F_a(i, j, k) = 3.0 * (F_a(i, j-1, k) - F_a(i, j-2, k)) + F_a(i, j-3, k);

          }
        }
      }
    }
    if(l_idx == 4)
    {
      for(int k = upper[2]; k >= lower[2]; k--)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            F_a(i, j, k) = 3.0 * (F_a(i, j, k+1) - F_a(i, j, k+2)) + F_a(i, j, k+3);

          }
        }
      }
    }
    if(l_idx == 5)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            F_a(i, j, k) = 3.0 * (F_a(i, j, k-1) - F_a(i, j, k-2)) + F_a(i, j, k-3);

          }
        }
      }
    }

    
  }

//   /************************updating codim = 2 boundaries****************/
//   codim = 2;

//   const std::vector<hier::BoundaryBox> & codim2_boxes =
//     geom->getCodimensionBoundaries(codim);

//   const idx_t n_codim2_boxes = static_cast<idx_t>(codim2_boxes.size());




//   for(int l = 0 ; l < n_codim2_boxes; l++)
//   {
//     hier::Box  boundary_fill_box =
//       geom->getBoundaryFillBox(
//         codim2_boxes[l], patch_box, F_a_pdata->getGhostCellWidth());
    
//     if(boundary_fill_box.empty()) continue;  


//     idx_t l_idx = codim2_boxes[l].getLocationIndex();
    
//     BSSNData bd = {0};

//     std::vector<idx_t> shift_vec;
//     if(l_idx == 0 || l_idx == 2 || l_idx == 4 || l_idx == 6)
//       shift_vec.push_back(STENCIL_ORDER_WIDTH);
//     else if(l_idx == 1 || l_idx == 3 || l_idx == 5 || l_idx == 7)
//       shift_vec.push_back(-STENCIL_ORDER_WIDTH);
//     else
//       shift_vec.push_back(0);
    
//     if(l_idx == 0 || l_idx == 1 || l_idx == 8 || l_idx == 10)
//       shift_vec.push_back(STENCIL_ORDER_WIDTH);
//     else if(l_idx == 2 || l_idx == 3 || l_idx ==9 || l_idx == 11)
//      shift_vec.push_back(-STENCIL_ORDER_WIDTH);
//     else
//       shift_vec.push_back(0);

//     if( l_idx == 4 || l_idx == 5 || l_idx == 8 || l_idx == 9)
//       shift_vec.push_back(STENCIL_ORDER_WIDTH);
//     else if(l_idx == 6 || l_idx == 7 || l_idx == 10 || l_idx == 11)
//       shift_vec.push_back(-STENCIL_ORDER_WIDTH);
//     else
//       shift_vec.push_back(0);

//     boundary_fill_box.shift(hier::IntVector(shift_vec));

//     boundary_fill_box *= patch_box;

//     //initialize dx for each patch
//     const real_t * dx = &(patch_geom->getDx())[0];

//     const idx_t * lower = &boundary_fill_box.lower()[0];
//     const idx_t * upper = &boundary_fill_box.upper()[0];

    
//     for(int k = lower[2]; k <= upper[2]; k++)
//     {
//       for(int j = lower[1]; j <= upper[1]; j++)
//       {
//         for(int i = lower[0]; i <= upper[0]; i++)
//         {
//           real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0;
//           real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0;
//           real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0;
        
//           F_a(i, j, k) = sqrt(pw2(x - origin[0]) + pw2(y - origin[1]) + pw2(z - origin[2]))
//             - getFonBD(patch, x, y, z, g_box);
//         }
//       }
//     }
//   }

//   /************************updating codim = 3 boundaries****************/
//   codim = 3;

//   const std::vector<hier::BoundaryBox> & codim3_boxes =
//     geom->getCodimensionBoundaries(codim);

//   const idx_t n_codim3_boxes = static_cast<idx_t>(codim3_boxes.size());




//   for(int l = 0 ; l < n_codim3_boxes; l++)
//   {
//     hier::Box boundary_fill_box =
//       geom->getBoundaryFillBox(
//         codim3_boxes[l], patch_box, F_a_pdata->getGhostCellWidth());

//     if(boundary_fill_box.empty()) continue;

//     idx_t l_idx = codim3_boxes[l].getLocationIndex();

//     std::vector<idx_t> shift_vec;
    
//     if(l_idx == 0 || l_idx == 2 || l_idx == 4 || l_idx == 6)
//       shift_vec.push_back(STENCIL_ORDER_WIDTH);
//     else
//       shift_vec.push_back(-STENCIL_ORDER_WIDTH);
    
//     if(l_idx == 0 || l_idx == 1 || l_idx == 4 || l_idx == 5)
//       shift_vec.push_back(STENCIL_ORDER_WIDTH);
//     else
//       shift_vec.push_back(-STENCIL_ORDER_WIDTH);

//     if( l_idx == 0 || l_idx == 1 || l_idx == 2 || l_idx == 3)
//       shift_vec.push_back(STENCIL_ORDER_WIDTH);
//     else
//       shift_vec.push_back(-STENCIL_ORDER_WIDTH);

//     boundary_fill_box.shift(hier::IntVector(shift_vec));

//     boundary_fill_box *= patch_box;

    
//     const idx_t * lower = &boundary_fill_box.lower()[0];
//     const idx_t * upper = &boundary_fill_box.upper()[0];

    
//     BSSNData bd = {0};


//     //initialize dx for each patch
//     const real_t * dx = &(patch_geom->getDx())[0];
  
//     for(int k = lower[2]; k <= upper[2]; k++)
//     {
//       for(int j = lower[1]; j <= upper[1]; j++)
//       {
//         for(int i = lower[0]; i <= upper[0]; i++)
//         {
//           real_t x = domain_lower[0] + (double)i * dx[0] + dx[0]/2.0;
//           real_t y = domain_lower[1] + (double)j * dx[1] + dx[1]/2.0;
//           real_t z = domain_lower[2] + (double)k * dx[2] + dx[2]/2.0;
        
//           F_a(i, j, k) = sqrt(pw2(x - origin[0]) + pw2(y - origin[1]) + pw2(z - origin[2]))
//             - getFonBD(patch, x, y, z, g_box);
//         }
//       }
//     }
//   }

// }
}
  
}
