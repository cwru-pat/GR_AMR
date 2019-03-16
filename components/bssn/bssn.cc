#include "bssn.h"
#include "../../cosmo_includes.h"
#include "../../utils/math.h"


using namespace SAMRAI;

namespace cosmo
{

/**
 * @brief Constructor for BSSN class
 */
BSSN::BSSN(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  const tbox::Dimension& dim_in,
  std::shared_ptr<tbox::Database> database_in,
  std::ostream* l_stream_in,
  real_t KO_damping_coefficient_in):
  lstream(l_stream_in),
  cosmo_bssn_db(database_in),
  dim(dim_in),
  KO_damping_coefficient(KO_damping_coefficient_in),
  gaugeHandler(new BSSNGaugeHandler(cosmo_bssn_db)),
  gd_eta(cosmo_bssn_db->getDoubleWithDefault("gd_eta",1.0)),
  normalize_Aij(cosmo_bssn_db->getBoolWithDefault("normalize_Aij",false)),
  normalize_gammaij(cosmo_bssn_db->getBoolWithDefault("normalize_gammaij",false)),
  Z4c_K1_DAMPING_AMPLITUDE(cosmo_bssn_db->getDoubleWithDefault("z4c_k1", 0.0)),
  Z4c_K2_DAMPING_AMPLITUDE(cosmo_bssn_db->getDoubleWithDefault("z4c_k2", 0.0)),
  chi_lower_bd(cosmo_bssn_db->getDoubleWithDefault("chi_lower_bd", 0)),
  alpha_lower_bd_for_L2(cosmo_bssn_db->getDoubleWithDefault("alpha_lower_bd_for_L2", 0.3)),
  K0(cosmo_bssn_db->getDoubleWithDefault("K0", 0))
{
  if(!USE_Z4C)
    Z4c_K1_DAMPING_AMPLITUDE = Z4c_K2_DAMPING_AMPLITUDE = 0;
  
  BSSN_APPLY_TO_FIELDS(VAR_INIT);
  BSSN_APPLY_TO_SOURCES(VAR_INIT);
  BSSN_APPLY_TO_GEN1_EXTRAS(VAR_INIT);

  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

  // creating variable context
  // (if corresponding context doe not exit, it creats one automatically)
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
#if USE_BACKUP_FIELDS
  std::shared_ptr<hier::VariableContext> context_b(
    variable_db->getContext("BACKUP"));
#endif

  // creating BSSN fields with contexts 
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_scratch, s, STENCIL_ORDER);
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_previous, p, STENCIL_ORDER);
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_active, a, STENCIL_ORDER);
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k1, k1, STENCIL_ORDER);
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k2, k2, STENCIL_ORDER);
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k3, k3, STENCIL_ORDER);
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k4, k4, STENCIL_ORDER);
#if USE_BACKUP_FIELDS
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_b, b, STENCIL_ORDER);
#endif

  // creating source fields only with ACTIVE context
  BSSN_APPLY_TO_SOURCES_ARGS(REG_TO_CONTEXT, context_active, a, STENCIL_ORDER);

  // creating extra fields with with ACTIVE context
  BSSN_APPLY_TO_GEN1_EXTRAS_ARGS(REG_TO_CONTEXT, context_active, a, STENCIL_ORDER);

  init(hierarchy);  
}

BSSN::~BSSN()
{
  
}

void BSSN::set_norm(
  const std::shared_ptr<hier::Patch>& patch, bool need_init_arr)
{

}
/**
 * @brief  normalize Aij or gammaij or both
 */
void BSSN::set_norm(
  const std::shared_ptr<hier::PatchLevel>& level)
{
  if(normalize_Aij == 0 && normalize_gammaij == 0) return;
  for (hier::PatchLevel::iterator p(level->begin());
       p != level->end(); ++p)
  {
    const std::shared_ptr<hier::Patch>& patch = *p;

    set_norm(patch, true);
  }
    /*
   * On all but the finest level, assign 0 to vector
   * weight to cells covered by finer cells.
   */

}



void BSSN::rescale_lapse(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int weight_idx)
{
}
  
/**
 * @brief  setting length of physical domain and chi lower bound
 * 
 */
void BSSN::init(const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  const double * domain_lower = &grid_geometry.getXLower()[0];
  const double * domain_upper = &grid_geometry.getXUpper()[0];

  for(int i = 0 ; i < DIM; i++)
    L[i] = domain_upper[i] - domain_lower[i];

  const double * dx = &grid_geometry.getDx()[0];

  std::string chi_lower_bd_type =
    cosmo_bssn_db->getStringWithDefault("chi_lower_bd_type", "");
  if(chi_lower_bd_type == "static_blackhole")
    chi_lower_bd = pow((dx[0] / (1<<(hierarchy->getMaxNumberOfLevels()-1)))/4.0, 4.0) / 10.0;
  
}

void BSSN::allocField(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln)
{
  std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_ALLOC);
}

void BSSN::allocSrc(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln)
{
  std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
  BSSN_APPLY_TO_SOURCES(EXTRA_ARRAY_ALLOC);
}
  
void BSSN::allocGen1(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln) 
{
  std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
  BSSN_APPLY_TO_GEN1_EXTRAS(EXTRA_ARRAY_ALLOC);
}

void BSSN::clearField(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  for(idx_t ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
  {
    clearField(hierarchy, ln);
  }
}
  
void BSSN::clearField(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln)
{
  math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);
  BSSN_APPLY_TO_FIELDS_ARGS(RK4_ARRAY_ZERO, hcellmath);
}

  
void BSSN::clearSrc(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  for(idx_t ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
  {
    clearSrc(hierarchy, ln);
  }
}

void BSSN::clearSrc(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln)
{
    
  math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);
  BSSN_APPLY_TO_SOURCES_ARGS(EXTRA_ARRAY_ZERO, hcellmath);
}

void BSSN::clearGen1(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  for(idx_t ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
  {
    clearGen1(hierarchy, ln);
  }
}
  
void BSSN::clearGen1(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln) 
{
  math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);
  BSSN_APPLY_TO_GEN1_EXTRAS_ARGS(EXTRA_ARRAY_ZERO, hcellmath);
}

/**
 * @brief  adding all active components to a list 
 * 
 */
void BSSN::addFieldsToList(std::vector<idx_t> &list)
{
  BSSN_APPLY_TO_FIELDS_ARGS(ADD_TO_LIST, list);
}

void BSSN::setExtraFieldData()
{
  
}

/**
 * @brief  some initialization for each step, currently does nothing
 * 
 */
void BSSN::stepInit(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  setExtraFieldData(); // Set extra field information (eg. derived field data for gauge conditions)
  // if enabling backup fields _b, copy _p data to _b data as backup
  
}

/**
 * @brief  RK evolve physical boundary for particular patch
 * 
 */
void BSSN::RKEvolvePatchBD(
  const std::shared_ptr<hier::Patch> & patch,
  real_t dt)
{
  
  std::shared_ptr<hier::PatchGeometry> geom (patch->getPatchGeometry());

  idx_t codim = 1;

  // getting all codimension 1 boxes
  const std::vector<hier::BoundaryBox> & codim1_boxes =
    geom->getCodimensionBoundaries(codim);

  const idx_t n_codim1_boxes = static_cast<idx_t>(codim1_boxes.size());

  // if it has no codimention 1 boundary, it has no other type of boundaries
  if(n_codim1_boxes == 0) return;

  initPData(patch);
  initMDA(patch);

  const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom( 
    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

    
  const hier::Box& patch_box = patch->getBox();


  for(int l = 0 ; l < n_codim1_boxes; l++)
  {
    hier::Box boundary_fill_box =
      geom->getBoundaryFillBox(
        codim1_boxes[l], patch_box, DIFFchi_a_pdata->getGhostCellWidth());

    if(boundary_fill_box.empty()) continue;


    idx_t l_idx = codim1_boxes[l].getLocationIndex();
    

    boundary_fill_box.shift(
      (hier::Box::dir_t)l_idx/2,
      (l_idx%2)?(-STENCIL_ORDER_WIDTH):STENCIL_ORDER_WIDTH);


    boundary_fill_box *= patch_box;

    //initialize dx for each patch
    const real_t * dx = &(patch_geom->getDx())[0];

    const idx_t * lower = &boundary_fill_box.lower()[0];
    const idx_t * upper = &boundary_fill_box.upper()[0];

    #pragma omp parallel for collapse(2)
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          BSSNData bd = {0};

          set_bd_values_bd(i, j, k, &bd, dx);
          BSSN_RK_EVOLVE_BD;
        }
      }
    }
  }
  /************************updating codim = 2 boundaries****************/
  codim = 2;

  const std::vector<hier::BoundaryBox> & codim2_boxes =
    geom->getCodimensionBoundaries(codim);

  const idx_t n_codim2_boxes = static_cast<idx_t>(codim2_boxes.size());




  for(int l = 0 ; l < n_codim2_boxes; l++)
  {
    hier::Box  boundary_fill_box =
      geom->getBoundaryFillBox(
        codim2_boxes[l], patch_box, DIFFchi_a_pdata->getGhostCellWidth());
    
    if(boundary_fill_box.empty()) continue;  


    idx_t l_idx = codim2_boxes[l].getLocationIndex();
    

    std::vector<idx_t> shift_vec;
    if(l_idx == 0 || l_idx == 2 || l_idx == 4 || l_idx == 6)
      shift_vec.push_back(STENCIL_ORDER_WIDTH);
    else if(l_idx == 1 || l_idx == 3 || l_idx == 5 || l_idx == 7)
      shift_vec.push_back(-STENCIL_ORDER_WIDTH);
    else
      shift_vec.push_back(0);
    
    if(l_idx == 0 || l_idx == 1 || l_idx == 8 || l_idx == 10)
      shift_vec.push_back(STENCIL_ORDER_WIDTH);
    else if(l_idx == 2 || l_idx == 3 || l_idx ==9 || l_idx == 11)
     shift_vec.push_back(-STENCIL_ORDER_WIDTH);
    else
      shift_vec.push_back(0);

    if( l_idx == 4 || l_idx == 5 || l_idx == 8 || l_idx == 9)
      shift_vec.push_back(STENCIL_ORDER_WIDTH);
    else if(l_idx == 6 || l_idx == 7 || l_idx == 10 || l_idx == 11)
      shift_vec.push_back(-STENCIL_ORDER_WIDTH);
    else
      shift_vec.push_back(0);

    boundary_fill_box.shift(hier::IntVector(shift_vec));

    boundary_fill_box *= patch_box;

    //initialize dx for each patch
    const real_t * dx = &(patch_geom->getDx())[0];

    const idx_t * lower = &boundary_fill_box.lower()[0];
    const idx_t * upper = &boundary_fill_box.upper()[0];

#pragma omp parallel for collapse(2)    
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          BSSNData bd = {0};
          set_bd_values_bd(i, j, k, &bd, dx);
          BSSN_RK_EVOLVE_BD;
        }
      }
    }
  }

  /************************updating codim = 3 boundaries****************/
  codim = 3;

  const std::vector<hier::BoundaryBox> & codim3_boxes =
    geom->getCodimensionBoundaries(codim);

  const idx_t n_codim3_boxes = static_cast<idx_t>(codim3_boxes.size());




  for(int l = 0 ; l < n_codim3_boxes; l++)
  {
    hier::Box boundary_fill_box =
      geom->getBoundaryFillBox(
        codim3_boxes[l], patch_box, DIFFchi_a_pdata->getGhostCellWidth());

    if(boundary_fill_box.empty()) continue;

    idx_t l_idx = codim3_boxes[l].getLocationIndex();

    std::vector<idx_t> shift_vec;
    
    if(l_idx == 0 || l_idx == 2 || l_idx == 4 || l_idx == 6)
      shift_vec.push_back(STENCIL_ORDER_WIDTH);
    else
      shift_vec.push_back(-STENCIL_ORDER_WIDTH);
    
    if(l_idx == 0 || l_idx == 1 || l_idx == 4 || l_idx == 5)
      shift_vec.push_back(STENCIL_ORDER_WIDTH);
    else
      shift_vec.push_back(-STENCIL_ORDER_WIDTH);

    if( l_idx == 0 || l_idx == 1 || l_idx == 2 || l_idx == 3)
      shift_vec.push_back(STENCIL_ORDER_WIDTH);
    else
      shift_vec.push_back(-STENCIL_ORDER_WIDTH);

    boundary_fill_box.shift(hier::IntVector(shift_vec));

    boundary_fill_box *= patch_box;

    const idx_t * lower = &boundary_fill_box.lower()[0];
    const idx_t * upper = &boundary_fill_box.upper()[0];

    


    //initialize dx for each patch
    const real_t * dx = &(patch_geom->getDx())[0];
    #pragma omp parallel for collapse(2)
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          BSSNData bd = {0};
          set_bd_values_bd(i, j, k, &bd, dx);
          BSSN_RK_EVOLVE_BD;
        }
      }
    }
  }
  
}

/**
 * @brief RK evolve patch interior
 */
void BSSN::RKEvolvePatch(
  const std::shared_ptr<hier::Patch> & patch, real_t dt)
{
  // might not need this function
  initPData(patch);
  initMDA(patch);

  const hier::Box& box = patch->getBox();
  
  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
  const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(  
    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  //initialize dx for each patch
  const real_t * dx = &(patch_geom->getDx())[0];


  bool flag = 0;
#pragma omp parallel for collapse(2)
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        BSSNData bd = {0};
        set_bd_values(i, j, k, &bd, dx);
        BSSN_RK_EVOLVE_PT;
      }
    }
  }

  return;
}

void BSSN::deBug(  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);
  const hier::Box& box = patch->getBox();
  
  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  BSSNData bd = {0};
  
  const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(  
    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  //initialize dx for each patch
  const real_t * dx = &(patch_geom->getDx())[0];


  
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        BSSN_APPLY_TO_FIELDS(BSSN_DEBUG);
      }
    }
  }
  return;
}

/**
 * @brief evolve fields on one cell
 * 
 */
void BSSN::RKEvolvePt(
  idx_t i, idx_t j, idx_t k, BSSNData &bd, const real_t dx[], real_t dt)
{
  set_bd_values(i, j, k, &bd, dx);
  BSSN_RK_EVOLVE_PT;
}

void BSSN::RKEvolvePtBd(
  idx_t i, idx_t j, idx_t k, BSSNData &bd, const real_t dx[], real_t dt,
  int l_idx, int codim)
{
  set_bd_values_bd(i, j, k, &bd, dx);
  BSSN_RK_EVOLVE_BD;
}


/**
 * @brief calculate K1 value on coarser level
 * 
 */
void BSSN::prepareForK1(
  const std::shared_ptr<hier::PatchLevel> & level,
  real_t to_t)
{
  if(level == NULL) return;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
    
    initPData(patch);
    initMDA(patch);

    const hier::Box& box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];


    // right branch of the tree
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(DIFFchi_a_idx)->getTime()) < EPS)
    {
#pragma omp parallel for collapse(2)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(RK4_INIT_R_K1);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(DIFFchi_p_idx)->getTime())
       - (patch->getPatchData(DIFFchi_a_idx)->getTime()
          - patch->getPatchData(DIFFchi_p_idx)->getTime())) < EPS)
    {
#pragma omp parallel for collapse(2)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(RK4_INIT_L_K1);
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
void BSSN::prepareForK2(
  const std::shared_ptr<hier::PatchLevel> & level,
  real_t to_t)
{
  if(level == NULL) return;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
    initPData(patch);
    initMDA(patch);

    const hier::Box& box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];


    // right branch of the tree
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(DIFFchi_a_idx)->getTime()) < EPS)
    {
#pragma omp parallel for collapse(2)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(RK4_INIT_R_K2);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(DIFFchi_p_idx)->getTime())
       - (patch->getPatchData(DIFFchi_a_idx)->getTime()
          - patch->getPatchData(DIFFchi_p_idx)->getTime())) < EPS)
    {
#pragma omp parallel for collapse(2)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(RK4_INIT_L_K2);
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
void BSSN::prepareForK3(
  const std::shared_ptr<hier::PatchLevel> & level,
  real_t to_t)
{
  if(level == NULL) return;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
    initPData(patch);
    initMDA(patch);

    const hier::Box& box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    // right branch of the tree
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(DIFFchi_a_idx)->getTime()) < EPS)
    {
        #pragma omp parallel for collapse(2)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(RK4_INIT_R_K3);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(DIFFchi_p_idx)->getTime())
       - (patch->getPatchData(DIFFchi_a_idx)->getTime()
          - patch->getPatchData(DIFFchi_p_idx)->getTime())) < EPS)
    {
        #pragma omp parallel for collapse(2)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(RK4_INIT_L_K3);
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
void BSSN::prepareForK4(
  const std::shared_ptr<hier::PatchLevel> & level,
  double to_t)
{
  if(level == NULL) return;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
    initPData(patch);
    initMDA(patch);

    const hier::Box& box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    
    


    // right branch of the tree
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(DIFFchi_a_idx)->getTime()) < EPS)
    {
        #pragma omp parallel for collapse(2)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(RK4_INIT_R_K4);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(DIFFchi_p_idx)->getTime())
       - (patch->getPatchData(DIFFchi_a_idx)->getTime()
          - patch->getPatchData(DIFFchi_p_idx)->getTime())) < EPS)
    {
        #pragma omp parallel for collapse(2)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(RK4_INIT_L_K4);
          }
        }
      }
    }
  }
}  

/**
 * @brief register "active" component of the fields to refiner
 * 
 * @param refiner
 * @param refiner operator
 */
void BSSN::registerRKRefinerActive(
  xfer::RefineAlgorithm& refiner,
  std::shared_ptr<hier::RefineOperator> &space_refine_op)
{
  BSSN_APPLY_TO_FIELDS_ARGS(REGISTER_SPACE_REFINE_A, refiner, space_refine_op);
}

  
/**
 * @brief register "scratch" component of the fields to refiner
 * 
 * @param refiner
 * @param refiner operator
 */
void BSSN::registerRKRefiner(
  xfer::RefineAlgorithm& refiner,
  std::shared_ptr<hier::RefineOperator> &space_refine_op)
{
  BSSN_APPLY_TO_FIELDS_ARGS(REGISTER_SPACE_REFINE_S, refiner, space_refine_op);
}


/**
 * @brief register "active" component of the fields to coarsener 
 * 
 * @param coarsener
 * @param coarse operator
 */
void BSSN::registerCoarsenActive(
  xfer::CoarsenAlgorithm& coarsener,
  std::shared_ptr<hier::CoarsenOperator>& coarsen_op)
{
  BSSN_APPLY_TO_FIELDS_ARGS(REGISTER_COARSEN_A, coarsener, coarsen_op);
}





/**
 * @brief copy active component to previous component of all BSSN fields
 * 
 */
void BSSN::copyAToP(
  math::HierarchyCellDataOpsReal<real_t> & hcellmath)
{
  BSSN_APPLY_TO_FIELDS(COPY_A_TO_P);
}

#if USE_BACKUP_FIELDS

void BSSN::copyPToB(
  math::HierarchyCellDataOpsReal<real_t> & hcellmath)
{
  BSSN_APPLY_TO_FIELDS(COPY_P_TO_B);
}

void BSSN::copyBToP(
  math::HierarchyCellDataOpsReal<real_t> & hcellmath)
{
  BSSN_APPLY_TO_FIELDS(COPY_B_TO_P);
}

void BSSN::copyBToA(
  math::HierarchyCellDataOpsReal<real_t> & hcellmath)
{
  BSSN_APPLY_TO_FIELDS(COPY_B_TO_A);
}


#endif


/**
 * @brief initilizing all pointers for every components's patchdata
 * 
 */
void BSSN::initPData(
  const std::shared_ptr<hier::Patch> & patch)
{
  BSSN_APPLY_TO_FIELDS(PDATA_ALL_INIT);
  BSSN_APPLY_TO_SOURCES_ARGS(PDATA_INIT, a);
  BSSN_APPLY_TO_GEN1_EXTRAS_ARGS(PDATA_INIT, a);
}


/**
 * @brief initilizing all pointers for every components's array access
 * 
 */
void BSSN::initMDA(
  const std::shared_ptr<hier::Patch> & patch)
{
  BSSN_APPLY_TO_FIELDS(MDA_ACCESS_ALL_INIT);

  BSSN_APPLY_TO_SOURCES_ARGS(MDA_ACCESS_INIT, a);
  
  BSSN_APPLY_TO_GEN1_EXTRAS_ARGS(MDA_ACCESS_INIT, a);
}

/**
 * @brief set the time of previous components of all BSSN fields as from_t
 *        set the time of active components of all BSSN fields as to_t
 */ 
void BSSN::setLevelTime(
  const std::shared_ptr<hier::PatchLevel> & level,
  double from_t, double to_t)
{
  BSSN_APPLY_TO_FIELDS_ARGS(SET_LEVEL_TIME, from_t, to_t);
}

/**
 * @brief  finalize k1 step for RK on both interior and boundary
 * 
 */
void BSSN::K1FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);
  
  const hier::Box& box = DIFFchi_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

    #pragma omp parallel for collapse(2)
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        BSSN_FINALIZE_K(1);
      }
    }
  }
}

void BSSN::K2FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  const hier::Box& box = DIFFchi_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  #pragma omp parallel for collapse(2)  
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        BSSN_FINALIZE_K(2);
      }
    }
  }
}


void BSSN::K3FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  const hier::Box& box = DIFFchi_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

    #pragma omp parallel for collapse(2)
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        BSSN_FINALIZE_K(3);
      }
    }
  }
}

// Only use this when you need to take the changes
// of some field as stop critieria
void BSSN::K4FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch, int ln, int max_ln)
{
  initPData(patch);
  initMDA(patch);

  //Only evolv inner grids for K4
  const hier::Box& box = DIFFchi_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  const double *dx = &patch_geom->getDx()[0];
  
  #pragma omp parallel for collapse(2)  
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        BSSN_FINALIZE_K(4);
      }
    }
  }
}


void BSSN::K4FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  //Only evolv inner grids for K4
  const hier::Box& box = DIFFchi_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  #pragma omp parallel for collapse(2)  
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        BSSN_FINALIZE_K(4);
      }
    }
  }
}



/*
******************************************************************************

 Calculate independent quantities for later use
(minimize # times derivatives are computed, etc)

******************************************************************************
*/

void BSSN::set_bd_values_bd(
  idx_t i, idx_t j, idx_t k, BSSNData *bd, const real_t dx[])
{
}

/**
 * @brief Populate values in a BSSNData struct
 * @details Compute all of them, except full metric m (TODO)
 * 
 * @param i x-index
 * @param j y-index
 * @param k z-index
 * @param bd BSSNData struct to populate
 */
void BSSN::set_bd_values(idx_t i, idx_t j, idx_t k, BSSNData *bd, const real_t dx[])
{
  bd->i = i;
  bd->j = j;
  bd->k = k;
  
  // need to set FRW quantities first

  //Have not figured out the initial value of FRW 
  
  bd->chi_FRW = 1;
  bd->K_FRW = 0;
  bd->rho_FRW = 0;
  bd->S_FRW = 0;

  bd->rho_P_avg = rho_P_avg;
  // draw data from cache
  set_local_vals(bd);
  bd->r        =   bd->DIFFr + bd->rho_FRW;
  bd->S        =   bd->DIFFS + bd->S_FRW;
  
  // pre-compute re-used quantities
  // gammas & derivs first

  // #if USE_EXPANSION
  // calculate_dexpN(bd,dx);
  // #endif
  // Christoffels depend on metric & derivs.
}


/**
 * @brief Set "local values"; set BSSNData values corresponding to field
 * values at a point.
 *
 * @param      bd    BSSNData struct with idx set.
 */
void BSSN::set_local_vals(BSSNData *bd)
{
  BSSN_APPLY_TO_FIELDS(RK4_SET_LOCAL_VALUES);
  BSSN_APPLY_TO_GEN1_EXTRAS(GEN1_SET_LOCAL_VALUES);
  BSSN_APPLY_TO_SOURCES(GEN1_SET_LOCAL_VALUES);
}

/**
 * @brief Compute and store inverse conformal difference metric components given
 * the conformal difference metric in a BSSNData struct
 * @details Computed assuming \f$det(\bar{\gamma}_{ij}) = 1\f$
 * 
 * @param i x-index
 * @param j y-index
 * @param k z-index
 * @param bd BSSNData containing initialized conformal difference metric components
 */
void BSSN::set_gammai_values(idx_t i, idx_t j, idx_t k, BSSNData *bd)
{
}

/**
 * @brief Calculate contravariant version of conformal trace-free extrinsic
 * curvature, \f$\bar{A}^{ij}\f$.
 *
 * @param bd BSSNData struct with inverse metric, Aij already computed.
 */

void BSSN::calculate_Acont(BSSNData *bd, const real_t dx[])
{
  // A^ij is calculated from A_ij by raising wrt. the conformal metric
}

/**
 * @brief Compute partial derivatives of the conformal metric, store in a
 * BSSNData instance
 *
 * @param bd BSSNData struct reference
 */
void BSSN::calculate_dgamma(BSSNData *bd, const real_t dx[])
{

}

/**
 * @brief Compute second partial derivatives of the conformal metric, store in
 * a BSSNData instance
 *
 * @param bd BSSNData struct reference
 */
void BSSN::calculate_ddgamma(BSSNData *bd, const real_t dx[])
{
}

/**
 * @brief Compute partial derivatives of the lapse and conformal factor, store
 * in a BSSNData instance
 *
 * @param bd BSSNData struct reference
 */
void BSSN::calculate_dalpha_dchi(BSSNData *bd, const real_t dx[])
{
}

/**
 * @brief Compute partial derivatives of the trace of the extrinsic curvature,
 * store in a BSSNData instance
 *
 * @param bd BSSNData struct reference
 */
void BSSN::calculate_dK(BSSNData *bd, const real_t dx[])
{
}

#if USE_Z4C
void BSSN::calculate_dtheta(BSSNData *bd, const real_t dx[])
{
}
#endif

#if USE_BSSN_SHIFT
void BSSN::calculate_dbeta(BSSNData *bd, const real_t dx[])
{
}
#endif

// #if USE_EXPANSION
// void BSSN::calculate_dexpN(BSSNData *bd, const real_t dx[])
// {
//   bd->d1expN = upwind_derivative(bd->i, bd->j, bd->k, 1, expN_a, dx, bd->beta1);
//   bd->d2expN = upwind_derivative(bd->i, bd->j, bd->k, 2, expN_a, dx, bd->beta2);
//   bd->d3expN = upwind_derivative(bd->i, bd->j, bd->k, 3, expN_a, dx, bd->beta3);
// }
// #endif


/*
******************************************************************************

Compute "dependent" quantities (depend on previously calc'd vals)

******************************************************************************
*/

void BSSN::calculate_conformal_christoffels(BSSNData *bd, const real_t dx[])
{
}

void BSSN::calculateDDchi(BSSNData *bd, const real_t dx[])
{
}


void BSSN::calculateDDalphaTF(BSSNData *bd, const real_t dx[])
{
  // double covariant derivatives - use non-unitary metric - extra pieces that depend on phi!
}

/* Calculate trace-free ricci tensor components */
void BSSN::calculateRicciTF(BSSNData *bd, const real_t dx[])
{
  return;
}




/*
******************************************************************************

Evolution equation calculations

******************************************************************************
*/

real_t BSSN::ev_DIFFgamma11(BSSNData *bd, const real_t dx[]) {return 0;}
real_t BSSN::ev_DIFFgamma12(BSSNData *bd, const real_t dx[]) {return 0;}
real_t BSSN::ev_DIFFgamma13(BSSNData *bd, const real_t dx[]) {return 0;}
real_t BSSN::ev_DIFFgamma22(BSSNData *bd, const real_t dx[]) {return 0;}
real_t BSSN::ev_DIFFgamma23(BSSNData *bd, const real_t dx[]) {return 0;}
real_t BSSN::ev_DIFFgamma33(BSSNData *bd, const real_t dx[]) {return 0;}

real_t BSSN::ev_A11(BSSNData *bd, const real_t dx[]) {return 0;}
real_t BSSN::ev_A12(BSSNData *bd, const real_t dx[]) {return 0;}
real_t BSSN::ev_A13(BSSNData *bd, const real_t dx[]) {return 0;}
real_t BSSN::ev_A22(BSSNData *bd, const real_t dx[]) {return 0;}
real_t BSSN::ev_A23(BSSNData *bd, const real_t dx[]) {return 0;}
real_t BSSN::ev_A33(BSSNData *bd, const real_t dx[]) {return 0;}

real_t BSSN::ev_Gamma1(BSSNData *bd, const real_t dx[]) {return 0;}
real_t BSSN::ev_Gamma2(BSSNData *bd, const real_t dx[]) {return 0;}
real_t BSSN::ev_Gamma3(BSSNData *bd, const real_t dx[]) {return 0;}


real_t BSSN::ev_DIFFK(BSSNData *bd, const real_t dx[])
{ return 0;;
}

  
real_t BSSN::ev_DIFFchi(BSSNData *bd, const real_t dx[])
{return 0;;
}

real_t BSSN::ev_DIFFalpha(BSSNData *bd, const real_t dx[])
{return 0;;
}


real_t BSSN::ev_theta(BSSNData *bd, const real_t dx[])
{return 0;
}


#if USE_BSSN_SHIFT
real_t BSSN::ev_beta1(BSSNData *bd, const real_t dx[])
{
  return gaugeHandler->ev_shift1(bd)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 1, beta1_a, dx, bd->beta1)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 2, beta1_a, dx, bd->beta2)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 3, beta1_a, dx, bd->beta3)
    - KO_dissipation_Q(bd->i, bd->j, bd->k, beta1_a, dx, KO_damping_coefficient);
}

real_t BSSN::ev_beta2(BSSNData *bd, const real_t dx[])
{
  return gaugeHandler->ev_shift2(bd)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 1, beta2_a, dx, bd->beta1)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 2, beta2_a, dx, bd->beta2)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 3, beta2_a, dx, bd->beta3)
    - KO_dissipation_Q(bd->i, bd->j, bd->k, beta2_a, dx, KO_damping_coefficient);
}

real_t BSSN::ev_beta3(BSSNData *bd, const real_t dx[])
{
  return gaugeHandler->ev_shift3(bd)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 1, beta3_a, dx, bd->beta1)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 2, beta3_a, dx, bd->beta2)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 3, beta3_a, dx, bd->beta3)
    - KO_dissipation_Q(bd->i, bd->j, bd->k, beta3_a, dx, KO_damping_coefficient);
}
#endif

#if USE_EXPANSION
real_t BSSN::ev_expN(BSSNData *bd, const real_t dx[])
{
  return
      upwind_derivative(bd->i, bd->j, bd->k, 1, expN_a, dx, bd->beta1)
    + upwind_derivative(bd->i, bd->j, bd->k, 2, expN_a, dx, bd->beta2)
    + upwind_derivative(bd->i, bd->j, bd->k, 3, expN_a, dx, bd->beta3)
    -bd->alpha * bd->K/3.0
    - KO_dissipation_Q(bd->i, bd->j, bd->k, expN_a, dx, KO_damping_coefficient);
}
#endif


#if USE_GAMMA_DRIVER
real_t BSSN::ev_auxB1(BSSNData *bd, const real_t dx[])
{
  return 0.75*ev_Gamma1(bd,dx)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 1, auxB1_a, dx, bd->beta1)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 2, auxB1_a, dx, bd->beta2)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 3, auxB1_a, dx, bd->beta3)
    - gd_eta * bd->auxB1
    - KO_dissipation_Q(bd->i, bd->j, bd->k, auxB1_a, dx, KO_damping_coefficient);
}

real_t BSSN::ev_auxB2(BSSNData *bd, const real_t dx[])
{
  return 0.75*ev_Gamma2(bd,dx)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 1, auxB2_a, dx, bd->beta1)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 2, auxB2_a, dx, bd->beta2)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 3, auxB2_a, dx, bd->beta3)
    - gd_eta * bd->auxB2
    - KO_dissipation_Q(bd->i, bd->j, bd->k, auxB2_a, dx, KO_damping_coefficient);
}

real_t BSSN::ev_auxB3(BSSNData *bd, const real_t dx[])
{
  return 0.75*ev_Gamma3(bd,dx)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 1, auxB3_a, dx, bd->beta1)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 2, auxB3_a, dx, bd->beta2)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 3, auxB3_a, dx, bd->beta3)
    - gd_eta * bd->auxB3
    - KO_dissipation_Q(bd->i, bd->j, bd->k, auxB3_a, dx, KO_damping_coefficient);
}
#endif

#if USE_PROPER_TIME
real_t BSSN::ev_tau(BSSNData *bd, const real_t dx[])
{
  return bd->alpha;
}

#endif

real_t BSSN::ev_H(BSSNData *bd, const real_t dx[])
{
  return -4.0 * PI * bd->rho_P_avg;
}

real_t BSSN::ev_a(BSSNData *bd, const real_t dx[])
{
  return bd->a * bd->H;
}

/*
******************************************************************************
Evolution equation calculations at boundary
******************************************************************************
*/

real_t BSSN::ev_DIFFgamma11_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}
real_t BSSN::ev_DIFFgamma12_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}
real_t BSSN::ev_DIFFgamma13_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}
real_t BSSN::ev_DIFFgamma22_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}
real_t BSSN::ev_DIFFgamma23_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}
real_t BSSN::ev_DIFFgamma33_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}


real_t BSSN::ev_A11_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}
real_t BSSN::ev_A12_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}
real_t BSSN::ev_A13_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}
real_t BSSN::ev_A22_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}
real_t BSSN::ev_A23_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}
real_t BSSN::ev_A33_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}

real_t BSSN::ev_Gamma1_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}
real_t BSSN::ev_Gamma2_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}
real_t BSSN::ev_Gamma3_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}


real_t BSSN::ev_DIFFK_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}

real_t BSSN::ev_DIFFchi_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}
real_t BSSN::ev_a_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}
real_t BSSN::ev_H_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}


real_t BSSN::ev_DIFFalpha_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}




real_t BSSN::ev_theta_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{return 0;
}


#if USE_BSSN_SHIFT
real_t BSSN::ev_beta1_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, beta1_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, beta1_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, beta1_a, dx, l_idx, codim) * bd->z 
            + bd->beta1           );
}

real_t BSSN::ev_beta2_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, beta2_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, beta2_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, beta2_a, dx, l_idx, codim) * bd->z 
            + bd->beta2           );
}

real_t BSSN::ev_beta3_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, beta3_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, beta3_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, beta3_a, dx, l_idx, codim) * bd->z 
            + bd->beta3           );
}

#if USE_EXPANSION
real_t BSSN::ev_expN_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, expN_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, expN_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, expN_a, dx, l_idx, codim) * bd->z 
            + bd->expN           );
}
#endif

#endif

#if USE_GAMMA_DRIVER
real_t BSSN::ev_auxB1_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, auxB1_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, auxB1_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, auxB1_a, dx, l_idx, codim) * bd->z 
            + bd->auxB1           );
}

real_t BSSN::ev_auxB2_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, auxB2_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, auxB2_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, auxB2_a, dx, l_idx, codim) * bd->z 
            + bd->auxB2           );

}

real_t BSSN::ev_auxB3_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, auxB3_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, auxB3_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, auxB3_a, dx, l_idx, codim) * bd->z 
            + bd->auxB3           );
}
#endif

#if USE_PROPER_TIME
real_t BSSN::ev_tau_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return bd->alpha;
}
#endif

void BSSN::output_max_H_constaint(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t weight_idx)
{
  return;
}

void BSSN::output_L2_H_constaint(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t weight_idx,   CosmoPatchStrategy * cosmoPS, double exclude_radius)
{
  return;
}


void BSSN::output_L2_H_constaint(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t weight_idx,   CosmoPatchStrategy * cosmoPS)
{
  return;
}


/*
******************************************************************************

Constraint violtion calculations

******************************************************************************
*/


real_t BSSN::hamiltonianConstraintCalc(BSSNData *bd, const real_t dx[])
{
}

real_t BSSN::hamiltonianConstraintScale(BSSNData *bd, const real_t dx[])
{
}


real_t BSSN::momentumConstraintCalc(BSSNData *bd, const real_t dx[])
{
}

real_t BSSN::momentumConstraintScale(BSSNData *bd, const real_t dx[])
{

}

} // namespace cosmo
