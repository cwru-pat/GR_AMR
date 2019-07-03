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
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_scratch, s, GHOST_WIDTH);
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_previous, p, GHOST_WIDTH);
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_active, a, GHOST_WIDTH);
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k1, k1, GHOST_WIDTH);
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k2, k2, GHOST_WIDTH);
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k3, k3, GHOST_WIDTH);
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k4, k4, GHOST_WIDTH);
#if USE_BACKUP_FIELDS
  BSSN_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_b, b, GHOST_WIDTH);
#endif

  // creating source fields only with ACTIVE context
  BSSN_APPLY_TO_SOURCES_ARGS(REG_TO_CONTEXT, context_active, a, GHOST_WIDTH);

  // creating extra fields with with ACTIVE context
  BSSN_APPLY_TO_GEN1_EXTRAS_ARGS(REG_TO_CONTEXT, context_active, a, GHOST_WIDTH);

  init(hierarchy);  
}

BSSN::~BSSN()
{
  
}

void BSSN::set_norm(
  const std::shared_ptr<hier::Patch>& patch, bool need_init_arr)
{
  if(need_init_arr)
  {
    initPData(patch);
    initMDA(patch);
  }
  std::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

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
        // 1 - det(1 + DiffGamma)
        if(normalize_gammaij)
        {
          real_t one_minus_det_gamma = -1.0*(
            DIFFgamma11_a(i,j,k) + DIFFgamma22_a(i,j,k) + DIFFgamma33_a(i,j,k)
            - pw2(DIFFgamma12_a(i,j,k)) - pw2(DIFFgamma13_a(i,j,k)) - pw2(DIFFgamma23_a(i,j,k))
            + DIFFgamma11_a(i,j,k)*DIFFgamma22_a(i,j,k)
            + DIFFgamma11_a(i,j,k)*DIFFgamma33_a(i,j,k)
            + DIFFgamma22_a(i,j,k)*DIFFgamma33_a(i,j,k)
            - pw2(DIFFgamma23_a(i,j,k))*DIFFgamma11_a(i,j,k)
            - pw2(DIFFgamma13_a(i,j,k))*DIFFgamma22_a(i,j,k)
            - pw2(DIFFgamma12_a(i,j,k))*DIFFgamma33_a(i,j,k)
            + 2.0*DIFFgamma12_a(i,j,k)*DIFFgamma13_a(i,j,k)*DIFFgamma23_a(i,j,k)
            + DIFFgamma11_a(i,j,k)*DIFFgamma22_a(i,j,k)*DIFFgamma33_a(i,j,k)
          );

          // accurately compute 1 - det(g)^(1/3), without roundoff error
          // = -( det(g)^(1/3) - 1 )
          // = -( exp{log[det(g)^(1/3)]} - 1 )
          // = -( expm1{log[det(g)]/3} )
          // = -expm1{log1p[-one_minus_det_gamma]/3.0}
          real_t one_minus_det_gamma_thirdpow = -1.0*expm1(log1p(-1.0*one_minus_det_gamma)/3.0);

          // Perform the equivalent of re-scaling the conformal metric so det(gamma) = 1
          // gamma -> gamma / det(gamma)^(1/3)
          // DIFFgamma -> (delta + DiffGamma) / det(gamma)^(1/3) - delta
          //            = ( DiffGamma + delta*[1 - det(gamma)^(1/3)] ) / ( 1 - [1 - det(1 + DiffGamma)^1/3] )
          DIFFgamma11_a(i,j,k) = (DIFFgamma11_a(i,j,k) + one_minus_det_gamma_thirdpow) / (1.0 - one_minus_det_gamma_thirdpow);
          DIFFgamma22_a(i,j,k) = (DIFFgamma22_a(i,j,k) + one_minus_det_gamma_thirdpow) / (1.0 - one_minus_det_gamma_thirdpow);
          DIFFgamma33_a(i,j,k) = (DIFFgamma33_a(i,j,k) + one_minus_det_gamma_thirdpow) / (1.0 - one_minus_det_gamma_thirdpow);
          DIFFgamma12_a(i,j,k) = (DIFFgamma12_a(i,j,k)) / (1.0 - one_minus_det_gamma_thirdpow);
          DIFFgamma13_a(i,j,k) = (DIFFgamma13_a(i,j,k)) / (1.0 - one_minus_det_gamma_thirdpow);
          DIFFgamma23_a(i,j,k) = (DIFFgamma23_a(i,j,k)) / (1.0 - one_minus_det_gamma_thirdpow);
        }
        if(normalize_Aij)
        {
          // re-scale A_ij / ensure it is trace-free
          // need inverse gamma for finding Tr(A)
          real_t gammai11 = 1.0 + DIFFgamma22_a(i,j,k) + DIFFgamma33_a(i,j,k) - pw2(DIFFgamma23_a(i,j,k)) + DIFFgamma22_a(i,j,k)*DIFFgamma33_a(i,j,k);
          real_t gammai22 = 1.0 + DIFFgamma11_a(i,j,k) + DIFFgamma33_a(i,j,k) - pw2(DIFFgamma13_a(i,j,k)) + DIFFgamma11_a(i,j,k)*DIFFgamma33_a(i,j,k);
          real_t gammai33 = 1.0 + DIFFgamma11_a(i,j,k) + DIFFgamma22_a(i,j,k) - pw2(DIFFgamma12_a(i,j,k)) + DIFFgamma11_a(i,j,k)*DIFFgamma22_a(i,j,k);
          real_t gammai12 = DIFFgamma13_a(i,j,k)*DIFFgamma23_a(i,j,k) - DIFFgamma12_a(i,j,k)*(1.0 + DIFFgamma33_a(i,j,k));
          real_t gammai13 = DIFFgamma12_a(i,j,k)*DIFFgamma23_a(i,j,k) - DIFFgamma13_a(i,j,k)*(1.0 + DIFFgamma22_a(i,j,k));
          real_t gammai23 = DIFFgamma12_a(i,j,k)*DIFFgamma13_a(i,j,k) - DIFFgamma23_a(i,j,k)*(1.0 + DIFFgamma11_a(i,j,k));
          real_t trA = gammai11*A11_a(i,j,k) + gammai22*A22_a(i,j,k) + gammai33*A33_a(i,j,k)
            + 2.0*(gammai12*A12_a(i,j,k) + gammai13*A13_a(i,j,k) + gammai23*A23_a(i,j,k));
          // A_ij -> ( A_ij - 1/3 gamma_ij A )
          A11_a(i,j,k) = ( A11_a(i,j,k) - 1.0/3.0*(1.0 + DIFFgamma11_a(i,j,k))*trA ) ;
          A22_a(i,j,k) = ( A22_a(i,j,k) - 1.0/3.0*(1.0 + DIFFgamma22_a(i,j,k))*trA ) ;
          A33_a(i,j,k) = ( A33_a(i,j,k) - 1.0/3.0*(1.0 + DIFFgamma33_a(i,j,k))*trA ) ;
          A12_a(i,j,k) = ( A12_a(i,j,k) - 1.0/3.0*DIFFgamma12_a(i,j,k)*trA ) ;
          A13_a(i,j,k) = ( A13_a(i,j,k) - 1.0/3.0*DIFFgamma13_a(i,j,k)*trA ) ;
          A23_a(i,j,k) = ( A23_a(i,j,k) - 1.0/3.0*DIFFgamma23_a(i,j,k)*trA ) ;
        }
      }
    }
  }

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

void BSSN::set_time_dependent_fields(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, double cur_t)
{
  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln ++)
  {
    std::shared_ptr <hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

    
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
    {
      const std::shared_ptr<hier::Patch> & patch = *pit;

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
            //DIFFgamma11_a(i, j, k) = DIFFgamma22_a(i, j, k) = DIFFgamma33_a(i, j, k)
            DIFFchi_a(i, j, k) = pow(1.5*cur_t + 1, -2.0/3.0) - 1.0;
          }
        }
      }

    }
  }

}


void BSSN::rescale_lapse(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int weight_idx)
{
  double max_alpha = -1.0;
  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln ++)
  {
    std::shared_ptr <hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

    
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
    {
      const std::shared_ptr<hier::Patch> & patch = *pit;

      initPData(patch);
      initMDA(patch);

      
      const hier::Box& box = DIFFchi_a_pdata->getGhostBox();
  
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];

      std::shared_ptr<pdat::CellData<double> > weight(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_idx)));
      

      arr_t weight_array =
        pdat::ArrayDataAccess::access<DIM, double>(
          weight->getArrayData());

      
#pragma omp parallel for collapse(2) reduction(max : max_alpha)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            if(weight_array(i, j, k) > 0)
              max_alpha = tbox::MathUtilities<double>::Max(max_alpha, DIFFalpha_a(i, j, k));
          }
        }
      }

    }
  }
  
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  mpi.Barrier();
  if (mpi.getSize() > 1) {
    mpi.AllReduce(&max_alpha, 1, MPI_MAX);
  }

  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln ++)
  {
    std::shared_ptr <hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

    
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
    {
      const std::shared_ptr<hier::Patch> & patch = *pit;

      initPData(patch);
      initMDA(patch);
      
      const hier::Box& box = DIFFchi_a_pdata->getBox();
  
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];

#pragma omp parallel for collapse(2)   
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
             DIFFalpha_a(i, j, k) = (DIFFalpha_a(i, j, k) + 1.0)/ (max_alpha + 1.0) - 1.0;
          }
        }
      }

    }
  }

  
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
#if USE_BACKUP_FIELDS
  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
  {
    math::HierarchyCellDataOpsReal<real_t> hcellmath(hierarchy,ln,ln);
    copyPToB(hcellmath);
  }
#endif
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
      (l_idx%2)?(-GHOST_WIDTH):GHOST_WIDTH);


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
      shift_vec.push_back(GHOST_WIDTH);
    else if(l_idx == 1 || l_idx == 3 || l_idx == 5 || l_idx == 7)
      shift_vec.push_back(-GHOST_WIDTH);
    else
      shift_vec.push_back(0);
    
    if(l_idx == 0 || l_idx == 1 || l_idx == 8 || l_idx == 10)
      shift_vec.push_back(GHOST_WIDTH);
    else if(l_idx == 2 || l_idx == 3 || l_idx ==9 || l_idx == 11)
     shift_vec.push_back(-GHOST_WIDTH);
    else
      shift_vec.push_back(0);

    if( l_idx == 4 || l_idx == 5 || l_idx == 8 || l_idx == 9)
      shift_vec.push_back(GHOST_WIDTH);
    else if(l_idx == 6 || l_idx == 7 || l_idx == 10 || l_idx == 11)
      shift_vec.push_back(-GHOST_WIDTH);
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
      shift_vec.push_back(GHOST_WIDTH);
    else
      shift_vec.push_back(-GHOST_WIDTH);
    
    if(l_idx == 0 || l_idx == 1 || l_idx == 4 || l_idx == 5)
      shift_vec.push_back(GHOST_WIDTH);
    else
      shift_vec.push_back(-GHOST_WIDTH);

    if( l_idx == 0 || l_idx == 1 || l_idx == 2 || l_idx == 3)
      shift_vec.push_back(GHOST_WIDTH);
    else
      shift_vec.push_back(-GHOST_WIDTH);

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
  bd->i = i;
  bd->j = j;
  bd->k = k;

  // need to set FRW quantities first
  bd->chi_FRW = 1;
  bd->K_FRW = 0;
  bd->rho_FRW = 0;
  bd->S_FRW = 0;
  
  // draw data from cache
  set_local_vals(bd);
  set_gammai_values(i, j, k, bd);

  // non-DIFF quantities
  bd->chi      =   bd->DIFFchi + bd->chi_FRW;
  bd->K        =   bd->DIFFK + bd->K_FRW;
  bd->gamma11  =   bd->DIFFgamma11 + 1.0;
  bd->gamma12  =   bd->DIFFgamma12;
  bd->gamma13  =   bd->DIFFgamma13;
  bd->gamma22  =   bd->DIFFgamma22 + 1.0;
  bd->gamma23  =   bd->DIFFgamma23;
  bd->gamma33  =   bd->DIFFgamma33 + 1.0;
  bd->r        =   bd->DIFFr + bd->rho_FRW;
  bd->S        =   bd->DIFFS + bd->S_FRW;
  bd->alpha    =   bd->DIFFalpha + 1.0;

  
  #if USE_SOMMERFIELD_BOUNDARY
  
  bd->x = (dx[0] * ((real_t)i + 0.5)) - L[0] / 2.0 ;
  bd->y = (dx[1] * ((real_t)j + 0.5)) - L[1] / 2.0 ;
  bd->z = (dx[2] * ((real_t)k + 0.5)) - L[2] / 2.0 ;
          
  bd->norm = sqrt(bd->x*bd->x + bd->y*bd->y + bd->z*bd->z);
  #endif
  
  if(bd->chi < chi_lower_bd) bd->chi = chi_lower_bd;
}

#if USE_COSMOTRACE
void BSSN::set_bd_values_for_ray_tracing(idx_t i, idx_t j, idx_t k, BSSNData *bd, const real_t dx[])
{
  for(int si = -2; si <= 2; si++)
    for(int sj = -2; sj <= 2; sj++)
      for(int sk = -2; sk <= 2; sk++)
      {
        hier::Index temp_idx(i+si, j+sj, k+sk);
        if(DIFFchi_a_pdata->getGhostBox().contains(temp_idx) == false)
          TBOX_ERROR("EEEEEEEEEEEEEEEEEEEE\n");
      }
  bd->i = i;
  bd->j = j;
  bd->k = k;

  set_local_vals(bd);

  // non-DIFF quantities
  bd->chi      =   bd->DIFFchi + 1.0;
  bd->K        =   bd->DIFFK + bd->K_FRW;
  bd->gamma11  =   bd->DIFFgamma11 + 1.0;
  bd->gamma12  =   bd->DIFFgamma12;
  bd->gamma13  =   bd->DIFFgamma13;
  bd->gamma22  =   bd->DIFFgamma22 + 1.0;
  bd->gamma23  =   bd->DIFFgamma23;
  bd->gamma33  =   bd->DIFFgamma33 + 1.0;
  bd->r        =   bd->DIFFr + bd->rho_FRW;
  bd->S        =   bd->DIFFS + bd->S_FRW;
  bd->alpha    =   bd->DIFFalpha + 1.0;

  calculate_dgamma(bd,dx);
  calculate_dalpha_dchi(bd,dx);
# if USE_BSSN_SHIFT
  calculate_dbeta(bd,dx);
# endif
  
}

#endif

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

  bd->K0 = K0;
  bd->K_avg = K_avg;
  
  // draw data from cache
  set_local_vals(bd);
  set_gammai_values(i, j, k, bd);

  // non-DIFF quantities
  bd->chi      =   bd->DIFFchi + bd->chi_FRW;
  bd->K        =   bd->DIFFK + bd->K_FRW;
  bd->gamma11  =   bd->DIFFgamma11 + 1.0;
  bd->gamma12  =   bd->DIFFgamma12;
  bd->gamma13  =   bd->DIFFgamma13;
  bd->gamma22  =   bd->DIFFgamma22 + 1.0;
  bd->gamma23  =   bd->DIFFgamma23;
  bd->gamma33  =   bd->DIFFgamma33 + 1.0;
  bd->r        =   bd->DIFFr + bd->rho_FRW;
  bd->S        =   bd->DIFFS + bd->S_FRW;
  bd->alpha    =   bd->DIFFalpha + 1.0;
  
  // pre-compute re-used quantities
  // gammas & derivs first
  calculate_Acont(bd,dx);
  calculate_dgamma(bd,dx);
  calculate_ddgamma(bd,dx);
  calculate_dalpha_dchi(bd,dx);
  calculate_dK(bd,dx);
# if USE_BSSN_SHIFT
  calculate_dbeta(bd,dx);
# endif

  // #if USE_EXPANSION
  // calculate_dexpN(bd,dx);
  // #endif
  // Christoffels depend on metric & derivs.
  calculate_conformal_christoffels(bd,dx);
  // DDw depend on christoffels, metric, and derivs
  calculateDDchi(bd,dx);
  calculateDDalphaTF(bd,dx);
  // Ricci depends on DDchi
  calculateRicciTF(bd,dx);

# if USE_Z4C
  calculate_dtheta(bd,dx);
# endif

  if(bd->chi < chi_lower_bd) bd->chi = chi_lower_bd;
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
  bd->gammai11 = 1.0 + bd->DIFFgamma22 + bd->DIFFgamma33 - pw2(bd->DIFFgamma23) + bd->DIFFgamma22*bd->DIFFgamma33;
  bd->gammai22 = 1.0 + bd->DIFFgamma11 + bd->DIFFgamma33 - pw2(bd->DIFFgamma13) + bd->DIFFgamma11*bd->DIFFgamma33;
  bd->gammai33 = 1.0 + bd->DIFFgamma11 + bd->DIFFgamma22 - pw2(bd->DIFFgamma12) + bd->DIFFgamma11*bd->DIFFgamma22;
  bd->gammai12 = bd->DIFFgamma13*bd->DIFFgamma23 - bd->DIFFgamma12*(1.0 + bd->DIFFgamma33);
  bd->gammai13 = bd->DIFFgamma12*bd->DIFFgamma23 - bd->DIFFgamma13*(1.0 + bd->DIFFgamma22);
  bd->gammai23 = bd->DIFFgamma12*bd->DIFFgamma13 - bd->DIFFgamma23*(1.0 + bd->DIFFgamma11);
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
  BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_ACONT)

  // calculate A_ij A^ij term
  AijAij_a(bd->i,bd->j,bd->k) = bd->Acont11*bd->A11 + bd->Acont22*bd->A22 + bd->Acont33*bd->A33
      + 2.0*(bd->Acont12*bd->A12 + bd->Acont13*bd->A13 + bd->Acont23*bd->A23);
  bd->AijAij = AijAij_a(bd->i,bd->j,bd->k);
}

/**
 * @brief Compute partial derivatives of the conformal metric, store in a
 * BSSNData instance
 *
 * @param bd BSSNData struct reference
 */
void BSSN::calculate_dgamma(BSSNData *bd, const real_t dx[])
{
  BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_DGAMMA)
}

/**
 * @brief Compute second partial derivatives of the conformal metric, store in
 * a BSSNData instance
 *
 * @param bd BSSNData struct reference
 */
void BSSN::calculate_ddgamma(BSSNData *bd, const real_t dx[])
{
  BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_DIDJGAMMA_PERMS)
}

/**
 * @brief Compute partial derivatives of the lapse and conformal factor, store
 * in a BSSNData instance
 *
 * @param bd BSSNData struct reference
 */
void BSSN::calculate_dalpha_dchi(BSSNData *bd, const real_t dx[])
{
  // normal derivatives of phi
  bd->d1chi = derivative(bd->i, bd->j, bd->k, 1, DIFFchi_a, dx);
  bd->d2chi = derivative(bd->i, bd->j, bd->k, 2, DIFFchi_a, dx);
  bd->d3chi = derivative(bd->i, bd->j, bd->k, 3, DIFFchi_a, dx);

  // second derivatives of phi
  bd->d1d1chi = double_derivative(bd->i, bd->j, bd->k, 1, 1, DIFFchi_a, dx);
  bd->d2d2chi = double_derivative(bd->i, bd->j, bd->k, 2, 2, DIFFchi_a, dx);
  bd->d3d3chi = double_derivative(bd->i, bd->j, bd->k, 3, 3, DIFFchi_a, dx);
  bd->d1d2chi = double_derivative(bd->i, bd->j, bd->k, 1, 2, DIFFchi_a, dx);
  bd->d1d3chi = double_derivative(bd->i, bd->j, bd->k, 1, 3, DIFFchi_a, dx);
  bd->d2d3chi = double_derivative(bd->i, bd->j, bd->k, 2, 3, DIFFchi_a, dx);
  
  // normal derivatives of alpha
  bd->d1a = derivative(bd->i, bd->j, bd->k, 1, DIFFalpha_a, dx);
  bd->d2a = derivative(bd->i, bd->j, bd->k, 2, DIFFalpha_a, dx);
  bd->d3a = derivative(bd->i, bd->j, bd->k, 3, DIFFalpha_a, dx);
}

/**
 * @brief Compute partial derivatives of the trace of the extrinsic curvature,
 * store in a BSSNData instance
 *
 * @param bd BSSNData struct reference
 */
void BSSN::calculate_dK(BSSNData *bd, const real_t dx[])
{
  // normal derivatives of K
  bd->d1K = derivative(bd->i, bd->j, bd->k, 1, DIFFK_a, dx);
  bd->d2K = derivative(bd->i, bd->j, bd->k, 2, DIFFK_a, dx);
  bd->d3K = derivative(bd->i, bd->j, bd->k, 3, DIFFK_a, dx);
}

#if USE_Z4C
void BSSN::calculate_dtheta(BSSNData *bd, const real_t dx[])
{
  // normal derivatives of phi
  bd->d1theta = derivative(bd->i, bd->j, bd->k, 1, theta_a, dx);
  bd->d2theta = derivative(bd->i, bd->j, bd->k, 2, theta_a, dx);
  bd->d3theta = derivative(bd->i, bd->j, bd->k, 3, theta_a, dx);
}
#endif

#if USE_BSSN_SHIFT
void BSSN::calculate_dbeta(BSSNData *bd, const real_t dx[])
{
  bd->d1beta1 = derivative(bd->i, bd->j, bd->k, 1, beta1_a, dx);
  bd->d1beta2 = derivative(bd->i, bd->j, bd->k, 1, beta2_a, dx);
  bd->d1beta3 = derivative(bd->i, bd->j, bd->k, 1, beta3_a, dx);
  bd->d2beta1 = derivative(bd->i, bd->j, bd->k, 2, beta1_a, dx);
  bd->d2beta2 = derivative(bd->i, bd->j, bd->k, 2, beta2_a, dx);
  bd->d2beta3 = derivative(bd->i, bd->j, bd->k, 2, beta3_a, dx);
  bd->d3beta1 = derivative(bd->i, bd->j, bd->k, 3, beta1_a, dx);
  bd->d3beta2 = derivative(bd->i, bd->j, bd->k, 3, beta2_a, dx);
  bd->d3beta3 = derivative(bd->i, bd->j, bd->k, 3, beta3_a, dx);
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
  // christoffel symbols: \Gamma^i_{jk} = Gijk
  BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_CHRISTOFFEL)
  // "lowered" christoffel symbols: \Gamma_{ijk} = GLijk
  BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_CHRISTOFFEL_LOWER)

  bd->Gammad1 = bd->G111*bd->gammai11 + bd->G122*bd->gammai22 + bd->G133*bd->gammai33
    + 2.0*(bd->G112*bd->gammai12 + bd->G113*bd->gammai13 + bd->G123*bd->gammai23);
  bd->Gammad2 = bd->G211*bd->gammai11 + bd->G222*bd->gammai22 + bd->G233*bd->gammai33
    + 2.0*(bd->G212*bd->gammai12 + bd->G213*bd->gammai13 + bd->G223*bd->gammai23);
  bd->Gammad3 = bd->G311*bd->gammai11 + bd->G322*bd->gammai22 + bd->G333*bd->gammai33
    + 2.0*(bd->G312*bd->gammai12 + bd->G313*bd->gammai13 + bd->G323*bd->gammai23);
}

void BSSN::calculateDDchi(BSSNData *bd, const real_t dx[])
{
  // double covariant derivatives, using unitary metric
  bd->D1D1chi = bd->d1d1chi - (bd->G111*bd->d1chi + bd->G211*bd->d2chi + bd->G311*bd->d3chi);
  bd->D2D2chi = bd->d2d2chi - (bd->G122*bd->d1chi + bd->G222*bd->d2chi + bd->G322*bd->d3chi);
  bd->D3D3chi = bd->d3d3chi - (bd->G133*bd->d1chi + bd->G233*bd->d2chi + bd->G333*bd->d3chi);

  bd->D1D2chi = bd->d1d2chi - (bd->G112*bd->d1chi + bd->G212*bd->d2chi + bd->G312*bd->d3chi);
  bd->D1D3chi = bd->d1d3chi - (bd->G113*bd->d1chi + bd->G213*bd->d2chi + bd->G313*bd->d3chi);
  bd->D2D3chi = bd->d2d3chi - (bd->G123*bd->d1chi + bd->G223*bd->d2chi + bd->G323*bd->d3chi);  
}


void BSSN::calculateDDalphaTF(BSSNData *bd, const real_t dx[])
{
  // double covariant derivatives - use non-unitary metric - extra pieces that depend on phi!
  // the gammaIldlphi are needed for the BSSN_CALCULATE_DIDJALPHA macro
  real_t gammai1ldlchi = bd->gammai11*bd->d1chi + bd->gammai12*bd->d2chi + bd->gammai13*bd->d3chi;
  real_t gammai2ldlchi = bd->gammai21*bd->d1chi + bd->gammai22*bd->d2chi + bd->gammai23*bd->d3chi;
  real_t gammai3ldlchi = bd->gammai31*bd->d1chi + bd->gammai32*bd->d2chi + bd->gammai33*bd->d3chi;
  // Calculates full (not trace-free) piece:
  BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_DIDJALPHA)

  // subtract trace (traced with full spatial metric but subtracted later)
  bd->DDaTR = bd->gammai11*bd->D1D1aTF + bd->gammai22*bd->D2D2aTF + bd->gammai33*bd->D3D3aTF
      + 2.0*(bd->gammai12*bd->D1D2aTF + bd->gammai13*bd->D1D3aTF + bd->gammai23*bd->D2D3aTF);
  bd->D1D1aTF -= (1.0/3.0)*bd->gamma11*bd->DDaTR;
  bd->D1D2aTF -= (1.0/3.0)*bd->gamma12*bd->DDaTR;
  bd->D1D3aTF -= (1.0/3.0)*bd->gamma13*bd->DDaTR;
  bd->D2D2aTF -= (1.0/3.0)*bd->gamma22*bd->DDaTR;
  bd->D2D3aTF -= (1.0/3.0)*bd->gamma23*bd->DDaTR;
  bd->D3D3aTF -= (1.0/3.0)*bd->gamma33*bd->DDaTR;

  // scale trace back (=> contracted with "real" metric)
  bd->DDaTR *= pw2(bd->chi);
}

/* Calculate trace-free ricci tensor components */
void BSSN::calculateRicciTF(BSSNData *bd, const real_t dx[])
{
  // unitary pieces
  BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_RICCI_UNITARY)

  /* calculate unitary Ricci scalar. */
  bd->unitRicci = bd->Uricci11*bd->gammai11 + bd->Uricci22*bd->gammai22 + bd->Uricci33*bd->gammai33
            + 2.0*(bd->Uricci12*bd->gammai12 + bd->Uricci13*bd->gammai13 + bd->Uricci23*bd->gammai23);

  /* Phi- contribution */
# if EXCLUDE_SECOND_ORDER_FRW
  TBOX_ERROR("EXCLUDE_SECOND_ORDER_FRW cannot be true");
  real_t expression = (
    bd->gammai11*bd->D1D1phi + bd->gammai22*bd->D2D2phi + bd->gammai33*bd->D3D3phi
    + 2.0*( bd->gammai12*bd->D1D2phi + bd->gammai13*bd->D1D3phi + bd->gammai23*bd->D2D3phi )
  );
  bd->ricci11 = bd->Uricci11 - 2.0*( bd->D1D1phi + bd->gamma11*(expression) );
  bd->ricci12 = bd->Uricci12 - 2.0*( bd->D1D2phi + bd->gamma12*(expression) );
  bd->ricci13 = bd->Uricci13 - 2.0*( bd->D1D3phi + bd->gamma13*(expression) );
  bd->ricci22 = bd->Uricci22 - 2.0*( bd->D2D2phi + bd->gamma22*(expression) );
  bd->ricci23 = bd->Uricci23 - 2.0*( bd->D2D3phi + bd->gamma23*(expression) );
  bd->ricci33 = bd->Uricci33 - 2.0*( bd->D3D3phi + bd->gamma33*(expression) );
# else
  real_t expression = (
    bd->gammai11*(bd->D1D1chi - 2.0*bd->d1chi*bd->d1chi/bd->chi)
    + bd->gammai22*(bd->D2D2chi - 2.0*bd->d2chi*bd->d2chi/bd->chi)
    + bd->gammai33*(bd->D3D3chi - 2.0*bd->d3chi*bd->d3chi/bd->chi)
    + 2.0*(
      bd->gammai12*(bd->D1D2chi - 2.0*bd->d1chi*bd->d2chi/bd->chi)
      + bd->gammai13*(bd->D1D3chi - 2.0*bd->d1chi*bd->d3chi/bd->chi)
      + bd->gammai23*(bd->D2D3chi - 2.0*bd->d2chi*bd->d3chi/bd->chi)
    )
  );

  bd->ricci11 = bd->Uricci11 + ( bd->D1D1chi + bd->gamma11*(expression) )/bd->chi;
  bd->ricci12 = bd->Uricci12 + ( bd->D1D2chi + bd->gamma12*(expression) )/bd->chi;
  bd->ricci13 = bd->Uricci13 + ( bd->D1D3chi + bd->gamma13*(expression) )/bd->chi;
  bd->ricci22 = bd->Uricci22 + ( bd->D2D2chi + bd->gamma22*(expression) )/bd->chi;
  bd->ricci23 = bd->Uricci23 + ( bd->D2D3chi + bd->gamma23*(expression) )/bd->chi;
  bd->ricci33 = bd->Uricci33 + ( bd->D3D3chi + bd->gamma33*(expression) )/bd->chi;
# endif

  /* calculate full Ricci scalar at this point */
  bd->ricci = pw2(bd->chi)*(
      bd->ricci11*bd->gammai11 + bd->ricci22*bd->gammai22 + bd->ricci33*bd->gammai33
      + 2.0*(bd->ricci12*bd->gammai12 + bd->ricci13*bd->gammai13 + bd->ricci23*bd->gammai23)
    );
  /* store ricci scalar here too. */
  ricci_a(bd->i,bd->j,bd->k) = bd->ricci;

  /* remove trace. Note that \bar{gamma}_{ij}*\bar{gamma}^{kl}R_{kl} = (unbarred gammas). */
  bd->ricciTF11 = bd->ricci11 - 1.0/(3.0*pw2(bd->chi))*bd->gamma11*bd->ricci;
  bd->ricciTF12 = bd->ricci12 - 1.0/(3.0*pw2(bd->chi))*bd->gamma12*bd->ricci;
  bd->ricciTF13 = bd->ricci13 - 1.0/(3.0*pw2(bd->chi))*bd->gamma13*bd->ricci;
  bd->ricciTF22 = bd->ricci22 - 1.0/(3.0*pw2(bd->chi))*bd->gamma22*bd->ricci;
  bd->ricciTF23 = bd->ricci23 - 1.0/(3.0*pw2(bd->chi))*bd->gamma23*bd->ricci;
  bd->ricciTF33 = bd->ricci33 - 1.0/(3.0*pw2(bd->chi))*bd->gamma33*bd->ricci;

  return;
}




/*
******************************************************************************

Evolution equation calculations

******************************************************************************
*/

real_t BSSN::ev_DIFFgamma11(BSSNData *bd, const real_t dx[]) { return BSSN_DT_DIFFGAMMAIJ(1, 1)  - KO_dissipation_Q(bd->i, bd->j, bd->k, DIFFgamma11_a, dx, KO_damping_coefficient); }
real_t BSSN::ev_DIFFgamma12(BSSNData *bd, const real_t dx[]) { return BSSN_DT_DIFFGAMMAIJ(1, 2)  - KO_dissipation_Q(bd->i, bd->j, bd->k, DIFFgamma12_a, dx, KO_damping_coefficient); }
real_t BSSN::ev_DIFFgamma13(BSSNData *bd, const real_t dx[]) { return BSSN_DT_DIFFGAMMAIJ(1, 3)  - KO_dissipation_Q(bd->i, bd->j, bd->k, DIFFgamma13_a, dx, KO_damping_coefficient); }
real_t BSSN::ev_DIFFgamma22(BSSNData *bd, const real_t dx[]) { return BSSN_DT_DIFFGAMMAIJ(2, 2)  - KO_dissipation_Q(bd->i, bd->j, bd->k, DIFFgamma22_a, dx, KO_damping_coefficient); }
real_t BSSN::ev_DIFFgamma23(BSSNData *bd, const real_t dx[]) { return BSSN_DT_DIFFGAMMAIJ(2, 3)  - KO_dissipation_Q(bd->i, bd->j, bd->k, DIFFgamma23_a, dx, KO_damping_coefficient); }
real_t BSSN::ev_DIFFgamma33(BSSNData *bd, const real_t dx[]) { return BSSN_DT_DIFFGAMMAIJ(3, 3)  - KO_dissipation_Q(bd->i, bd->j, bd->k, DIFFgamma33_a, dx, KO_damping_coefficient); }

real_t BSSN::ev_A11(BSSNData *bd, const real_t dx[]) { return BSSN_DT_AIJ(1, 1) - KO_dissipation_Q(bd->i, bd->j, bd->k, A11_a, dx, KO_damping_coefficient); }
real_t BSSN::ev_A12(BSSNData *bd, const real_t dx[]) { return BSSN_DT_AIJ(1, 2) - KO_dissipation_Q(bd->i, bd->j, bd->k, A12_a, dx, KO_damping_coefficient); }
real_t BSSN::ev_A13(BSSNData *bd, const real_t dx[]) { return BSSN_DT_AIJ(1, 3) - KO_dissipation_Q(bd->i, bd->j, bd->k, A13_a, dx, KO_damping_coefficient); }
real_t BSSN::ev_A22(BSSNData *bd, const real_t dx[]) { return BSSN_DT_AIJ(2, 2) - KO_dissipation_Q(bd->i, bd->j, bd->k, A22_a, dx, KO_damping_coefficient); }
real_t BSSN::ev_A23(BSSNData *bd, const real_t dx[]) { return BSSN_DT_AIJ(2, 3) - KO_dissipation_Q(bd->i, bd->j, bd->k, A23_a, dx, KO_damping_coefficient); }
real_t BSSN::ev_A33(BSSNData *bd, const real_t dx[]) { return BSSN_DT_AIJ(3, 3) - KO_dissipation_Q(bd->i, bd->j, bd->k, A33_a, dx, KO_damping_coefficient); }

real_t BSSN::ev_Gamma1(BSSNData *bd, const real_t dx[]) { return BSSN_DT_GAMMAI(1) - KO_dissipation_Q(bd->i, bd->j, bd->k, Gamma1_a, dx, KO_damping_coefficient); }
real_t BSSN::ev_Gamma2(BSSNData *bd, const real_t dx[]) { return BSSN_DT_GAMMAI(2) - KO_dissipation_Q(bd->i, bd->j, bd->k, Gamma2_a, dx, KO_damping_coefficient); }
real_t BSSN::ev_Gamma3(BSSNData *bd, const real_t dx[]) { return BSSN_DT_GAMMAI(3) - KO_dissipation_Q(bd->i, bd->j, bd->k, Gamma3_a, dx, KO_damping_coefficient); }


real_t BSSN::ev_DIFFK(BSSNData *bd, const real_t dx[])
{
  return (
    - bd->DDaTR
    + bd->alpha*(
#       if EXCLUDE_SECOND_ORDER_FRW
          1.0/3.0*(bd->DIFFK + 2.0*bd->theta)*2.0*bd->K_FRW
#       else
          pw2(bd->K + 2.0 * bd->theta)/3.0
#       endif

#       if !(EXCLUDE_SECOND_ORDER_SMALL)
          + bd->AijAij
#       endif
    )
    + 4.0*PI*bd->alpha*(bd->r + bd->S)
#if USE_BSSN_SHIFT
    + upwind_derivative(bd->i, bd->j, bd->k, 1, DIFFK_a, dx, bd->beta1)
    + upwind_derivative(bd->i, bd->j, bd->k, 2, DIFFK_a, dx, bd->beta2)
    + upwind_derivative(bd->i, bd->j, bd->k, 3, DIFFK_a, dx, bd->beta3)
#endif
    - bd->alpha*Z4c_K1_DAMPING_AMPLITUDE*(1.0 - Z4c_K2_DAMPING_AMPLITUDE)*bd->theta
    - KO_dissipation_Q(bd->i, bd->j, bd->k, DIFFK_a, dx, KO_damping_coefficient)
  );
}

  
real_t BSSN::ev_DIFFchi(BSSNData *bd, const real_t dx[])
{
  return (
    1.0/3.0*(
      bd->alpha* bd->chi * (bd->K + 2.0 * bd->theta)
      #if USE_BSSN_SHIFT
      - bd->chi * ( bd->d1beta1 + bd->d2beta2 + bd->d3beta3 )
      #endif
    )
    #if USE_BSSN_SHIFT
    + upwind_derivative(bd->i, bd->j, bd->k, 1, DIFFchi_a, dx, bd->beta1)
    + upwind_derivative(bd->i, bd->j, bd->k, 2, DIFFchi_a, dx, bd->beta2)
    + upwind_derivative(bd->i, bd->j, bd->k, 3, DIFFchi_a, dx, bd->beta3)
    #endif
    - KO_dissipation_Q(bd->i, bd->j, bd->k, DIFFchi_a, dx, KO_damping_coefficient)
  );
}

real_t BSSN::ev_DIFFalpha(BSSNData *bd, const real_t dx[])
{
  return gaugeHandler->ev_lapse(bd)
    #if USE_BSSN_SHIFT
    + upwind_derivative(bd->i, bd->j, bd->k, 1, DIFFalpha_a, dx, bd->beta1)
    + upwind_derivative(bd->i, bd->j, bd->k, 2, DIFFalpha_a, dx, bd->beta2)
    + upwind_derivative(bd->i, bd->j, bd->k, 3, DIFFalpha_a, dx, bd->beta3)
    #endif
    - KO_dissipation_Q(bd->i, bd->j, bd->k, DIFFalpha_a, dx, KO_damping_coefficient);
}


real_t BSSN::ev_theta(BSSNData *bd, const real_t dx[])
{
  #if USE_Z4C
  return (
    0.5*bd->alpha*(
      bd->ricci + 2.0/3.0*pw2(bd->K + 2.0 * bd->theta)
      - bd->AijAij - 16.0*PI* bd->r)
    #if USE_BSSN_SHIFT
    + upwind_derivative(bd->i, bd->j, bd->k, 1, theta_a, dx, bd->beta1)
    + upwind_derivative(bd->i, bd->j, bd->k, 2, theta_a, dx, bd->beta2)
    + upwind_derivative(bd->i, bd->j, bd->k, 3, theta_a, dx, bd->beta3)
    #endif
    - bd->alpha*Z4c_K1_DAMPING_AMPLITUDE*(2.0 + Z4c_K2_DAMPING_AMPLITUDE)*bd->theta
  ) - KO_dissipation_Q(bd->i, bd->j, bd->k, theta_a, dx, KO_damping_coefficient);
  #endif
  return 0;
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

/*
******************************************************************************
Evolution equation calculations at boundary
******************************************************************************
*/

real_t BSSN::ev_DIFFgamma11_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, DIFFgamma11_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, DIFFgamma11_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, DIFFgamma11_a, dx, l_idx, codim) * bd->z 
            + bd->DIFFgamma11           );
}
real_t BSSN::ev_DIFFgamma12_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, DIFFgamma12_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, DIFFgamma12_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, DIFFgamma12_a, dx, l_idx, codim) * bd->z 
            + bd->DIFFgamma12           );
}
real_t BSSN::ev_DIFFgamma13_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, DIFFgamma13_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, DIFFgamma13_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, DIFFgamma13_a, dx, l_idx, codim) * bd->z 
            + bd->DIFFgamma13           );
}
real_t BSSN::ev_DIFFgamma22_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, DIFFgamma22_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, DIFFgamma22_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, DIFFgamma22_a, dx, l_idx, codim) * bd->z 
            + bd->DIFFgamma22           );
}
real_t BSSN::ev_DIFFgamma23_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, DIFFgamma23_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, DIFFgamma23_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, DIFFgamma23_a, dx, l_idx, codim) * bd->z 
            + bd->DIFFgamma23           );
}
real_t BSSN::ev_DIFFgamma33_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, DIFFgamma33_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, DIFFgamma33_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, DIFFgamma33_a, dx, l_idx, codim) * bd->z 
            + bd->DIFFgamma33           );
}


real_t BSSN::ev_A11_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, A11_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, A11_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, A11_a, dx, l_idx, codim) * bd->z 
            + bd->A11           );
}
real_t BSSN::ev_A12_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, A12_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, A12_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, A12_a, dx, l_idx, codim) * bd->z 
            + bd->A12           );
}
real_t BSSN::ev_A13_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, A13_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, A13_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, A13_a, dx, l_idx, codim) * bd->z 
            + bd->A13           );
}
real_t BSSN::ev_A22_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, A22_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, A22_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, A22_a, dx, l_idx, codim) * bd->z 
            + bd->A22           );
}
real_t BSSN::ev_A23_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, A23_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, A23_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, A23_a, dx, l_idx, codim) * bd->z 
            + bd->A23           );
}
real_t BSSN::ev_A33_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, A33_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, A33_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, A33_a, dx, l_idx, codim) * bd->z 
            + bd->A33           );
}

real_t BSSN::ev_Gamma1_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, Gamma1_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, Gamma1_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, Gamma1_a, dx, l_idx, codim) * bd->z 
            + bd->Gamma1           );
}
real_t BSSN::ev_Gamma2_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, Gamma2_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, Gamma2_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, Gamma2_a, dx, l_idx, codim) * bd->z 
            + bd->Gamma2           );
}
real_t BSSN::ev_Gamma3_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, Gamma3_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, Gamma3_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, Gamma3_a, dx, l_idx, codim) * bd->z 
            + bd->Gamma3           );
}


real_t BSSN::ev_DIFFK_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, DIFFK_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, DIFFK_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, DIFFK_a, dx, l_idx, codim) * bd->z 
            + bd->DIFFK - bd->K0           );
}

real_t BSSN::ev_DIFFchi_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, DIFFchi_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, DIFFchi_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, DIFFchi_a, dx, l_idx, codim) * bd->z 
            + bd->DIFFchi           );
}
  
real_t BSSN::ev_DIFFalpha_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, DIFFalpha_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, DIFFalpha_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, DIFFalpha_a, dx, l_idx, codim) * bd->z 
            + bd->DIFFalpha           );
}




real_t BSSN::ev_theta_bd(BSSNData *bd, const real_t dx[], idx_t l_idx, idx_t codim)
{
#if USE_Z4C  
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, theta_a, dx, l_idx, codim) * bd->x
            + bd_derivative(bd->i, bd->j, bd->k, 2, theta_a, dx, l_idx, codim) * bd->y
            + bd_derivative(bd->i, bd->j, bd->k, 3, theta_a, dx, l_idx, codim) * bd->z 
            + bd->theta           );
#endif
  return 0;
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
  double max_H=0, max_H_scaled = 0;
  double max_M=0, max_M_scaled = 0;
  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln ++)
  {
    std::shared_ptr <hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
    {
      const std::shared_ptr<hier::Patch> & patch = *pit;

      const hier::Box& box = patch->getBox();

      const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      initPData(patch);

      initMDA(patch);


      std::shared_ptr<pdat::CellData<double> > weight(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_idx)));
      

      arr_t weight_array =
        pdat::ArrayDataAccess::access<DIM, double>(
          weight->getArrayData());


      
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];
      
      const double *dx = &patch_geom->getDx()[0];

      idx_t NX = round(L[0] / dx[0]);
      idx_t NY = round(L[1] / dx[1]); 
      idx_t NZ = round(L[2] / dx[2]);

      
#pragma omp parallel for collapse(2) reduction(max : max_H, max_M, max_H_scaled, max_M_scaled)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSNData bd = {0};

            set_bd_values(i,j,k,&bd,dx);
            double xx = (dx[0] * ((real_t)i + 0.5)) - L[0] / 2.0 ;
            double yy = (dx[1] * ((real_t)j + 0.5)) - L[1] / 2.0 ;
            double zz = (dx[2] * ((real_t)k + 0.5)) - L[2] / 2.0 ;
            double r = sqrt(pw2(xx) + pw2(yy) + pw2(zz));

            if(weight_array(i,j,k) > 0 )
            {
              max_H = tbox::MathUtilities<double>::Max(
                max_H, tbox::MathUtilities<double>::Abs(
                  hamiltonianConstraintCalc(&bd, dx)));
               max_H_scaled = tbox::MathUtilities<double>::Max(
                max_H_scaled, tbox::MathUtilities<double>::Abs(
                  hamiltonianConstraintCalc(&bd, dx)/
                  hamiltonianConstraintScale(&bd,dx)));

               max_M = tbox::MathUtilities<double>::Max(
                max_M, tbox::MathUtilities<double>::Abs(
                  momentumConstraintCalc(&bd, dx)));
               max_M_scaled = tbox::MathUtilities<double>::Max(
                max_M_scaled, tbox::MathUtilities<double>::Abs(
                  momentumConstraintCalc(&bd, dx)/
                  momentumConstraintScale(&bd,dx)));

            }
          }
        }
      }
     }
  }

  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  mpi.Barrier();
  if (mpi.getSize() > 1) {
    mpi.AllReduce(&max_H, 1, MPI_MAX);
    mpi.AllReduce(&max_H_scaled, 1, MPI_MAX);
    mpi.AllReduce(&max_M, 1, MPI_MAX);
    mpi.AllReduce(&max_M_scaled, 1, MPI_MAX);

  }

  tbox::pout<<"Max Hamiltonian constraint is "
            <<max_H<<"\n";
  tbox::pout<<"Max Momentum constraint is "
            <<max_M<<"\n";
  return;
}

void BSSN::output_L2_H_constaint(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t weight_idx,   CosmoPatchStrategy * cosmoPS, double exclude_radius)
{ 
  double H_L2=0;

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;
  const double * dx = &grid_geometry.getDx()[0];

  const int base_nx = round(L[0] / dx[0]);
  const int base_ny = round(L[1] / dx[1]);
  const int base_nz = round(L[2] / dx[2]);
  
  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln ++)
  {
    std::shared_ptr <hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

    
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
    {
      const std::shared_ptr<hier::Patch> & patch = *pit;

      const hier::Box& box = patch->getBox();

      const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      initPData(patch);

      initMDA(patch);


      std::shared_ptr<pdat::CellData<double> > weight(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_idx)));
      

      arr_t weight_array =
        pdat::ArrayDataAccess::access<DIM, double>(
          weight->getArrayData());


      
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];
      
      const double *dx = &patch_geom->getDx()[0];

      

#pragma omp parallel for collapse(2) reduction(+:H_L2)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            if(cosmoPS->is_time_dependent &&
               (i < GHOST_WIDTH || i >= base_nx - GHOST_WIDTH ||
                j < GHOST_WIDTH || j >= base_ny - GHOST_WIDTH ||
                k < GHOST_WIDTH || k >= base_nz - GHOST_WIDTH))
              continue;
            
            BSSNData bd = {0};

            set_bd_values(i,j,k,&bd,dx);

            double xx = (dx[0] * ((real_t)i + 0.5)) - L[0] / 2.0 ;
            double yy = (dx[1] * ((real_t)j + 0.5)) - L[1] / 2.0 ;
            double zz = (dx[2] * ((real_t)k + 0.5)) - L[2] / 2.0 ;
            double r = sqrt(pw2(xx) + pw2(yy) + pw2(zz));

            if(weight_array(i,j,k) > 0 && DIFFalpha_a(i,j,k) > alpha_lower_bd_for_L2 - 1.0
            && r > exclude_radius)
            {
              real_t h = hamiltonianConstraintCalc(&bd, dx);
              H_L2 += pw2(h) * weight_array(i,j,k) / (L[0] * L[1] * L[2]);
            }
          }
        }
      }
     }
  }
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  mpi.Barrier();
  if (mpi.getSize() > 1) {
    mpi.AllReduce(&H_L2, 1, MPI_SUM);
  }

  H_L2 = sqrt(H_L2);
  
  tbox::pout<<"L2 norm of Hamiltonian constraint is "<<H_L2<<"\n";

  
  return;
}


void BSSN::output_L2_H_constaint(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t weight_idx,   CosmoPatchStrategy * cosmoPS)
{
  double H_L2=0;

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;
  const double * dx = &grid_geometry.getDx()[0];

  const int base_nx = round(L[0] / dx[0]);
  const int base_ny = round(L[1] / dx[1]);
  const int base_nz = round(L[2] / dx[2]);
  
  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln ++)
  {
    std::shared_ptr <hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

    
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
    {
      const std::shared_ptr<hier::Patch> & patch = *pit;

      const hier::Box& box = patch->getBox();

      const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      initPData(patch);

      initMDA(patch);


      std::shared_ptr<pdat::CellData<double> > weight(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_idx)));
      

      arr_t weight_array =
        pdat::ArrayDataAccess::access<DIM, double>(
          weight->getArrayData());


      
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];
      
      const double *dx = &patch_geom->getDx()[0];

      

#pragma omp parallel for collapse(2) reduction(+:H_L2)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            if(cosmoPS->is_time_dependent &&
               (i < GHOST_WIDTH || i >= base_nx - GHOST_WIDTH ||
                j < GHOST_WIDTH || j >= base_ny - GHOST_WIDTH ||
                k < GHOST_WIDTH || k >= base_nz - GHOST_WIDTH))
              continue;
            
            BSSNData bd = {0};

            set_bd_values(i,j,k,&bd,dx);

            if(weight_array(i,j,k) > 0 && DIFFalpha_a(i,j,k) > alpha_lower_bd_for_L2 - 1.0)
            {
              real_t h = hamiltonianConstraintCalc(&bd, dx);
              H_L2 += pw2(h) * weight_array(i,j,k) / (L[0] * L[1] * L[2]);
            }
          }
        }
      }
     }
  }
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  mpi.Barrier();
  if (mpi.getSize() > 1) {
    mpi.AllReduce(&H_L2, 1, MPI_SUM);
  }

  H_L2 = sqrt(H_L2);
  
  tbox::pout<<"L2 norm of Hamiltonian constraint is "<<H_L2<<"\n";

  
  return;
}


/*
******************************************************************************

Constraint violtion calculations

******************************************************************************
*/


real_t BSSN::hamiltonianConstraintCalc(BSSNData *bd, const real_t dx[])
{
  return -pow(bd->chi, -2.5)/8.0*(
    bd->ricci + 2.0/3.0*pw2(bd->K + 2.0 * bd->theta) - bd->AijAij - 16.0*PI*bd->r
    );
}

real_t BSSN::hamiltonianConstraintScale(BSSNData *bd, const real_t dx[])
{

  // sqrt sum of sq. of terms for appx. mag / scale
  return 
    pow(bd->chi, -2.5)/8.0*sqrt( pw2(bd->ricci) + pw2(bd->AijAij)
          + pw2(2.0/3.0*pw2(bd->K + 2.0 * bd->theta))
          + pw2(16.0*PI*bd->r)
    );
}


real_t BSSN::momentumConstraintCalc(BSSNData *bd, const real_t dx[])
{
  real_t mi1 = BSSN_MI(1);
  real_t mi2 = BSSN_MI(2);
  real_t mi3 = BSSN_MI(3);


  return sqrt(pw2(mi1)+pw2(mi2)+pw2(mi3));
}

real_t BSSN::momentumConstraintScale(BSSNData *bd, const real_t dx[])
{
  // sqrt sum of sq. of terms for appx. mag / scale
  real_t mi1s = BSSN_MI_SCALE(1);
  real_t mi2s = BSSN_MI_SCALE(2);
  real_t mi3s = BSSN_MI_SCALE(3);

  return sqrt(pw2(mi1s)+pw2(mi2s)+pw2(mi3s));

}

} // namespace cosmo
