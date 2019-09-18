#include "dust_fluid.h"
#include "../../cosmo_includes.h"
#include "../../utils/math.h"


using namespace SAMRAI;

namespace cosmo
{

DustFluid::DustFluid(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  const tbox::Dimension& dim_in,
  std::shared_ptr<tbox::Database> database_in,
  std::ostream* l_stream_in):
  lstream(l_stream_in),
  cosmo_dust_fluid_db(database_in),
  dim(dim_in),
  is_test_fluid(cosmo_dust_fluid_db->getBoolWithDefault("is_test_fluid",false))
{
  DUST_FLUID_APPLY_TO_FIELDS(VAR_INIT);
  DUST_FLUID_APPLY_TO_DERIVED_FIELDS(VAR_INIT);

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

  
  DUST_FLUID_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_scratch, s, GHOST_WIDTH);
  DUST_FLUID_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_previous, p, GHOST_WIDTH);
  DUST_FLUID_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_active, a, GHOST_WIDTH);
  DUST_FLUID_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k1, k1, GHOST_WIDTH);
  DUST_FLUID_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k2, k2, GHOST_WIDTH);
  DUST_FLUID_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k3, k3, GHOST_WIDTH);
  DUST_FLUID_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k4, k4, GHOST_WIDTH);

  DUST_FLUID_APPLY_TO_DERIVED_FIELDS_ARGS(REG_TO_CONTEXT, context_active, a, GHOST_WIDTH);
}
DustFluid::~DustFluid()
{
  
}

void DustFluid::init(const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  return;
}
  
void DustFluid::alloc(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln)
{
  std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
  DUST_FLUID_APPLY_TO_FIELDS(RK4_ARRAY_ALLOC);
  DUST_FLUID_APPLY_TO_DERIVED_FIELDS(EXTRA_ARRAY_ALLOC);
  
}

void DustFluid::clear(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln)
{
  math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);
  DUST_FLUID_APPLY_TO_FIELDS_ARGS(RK4_ARRAY_ZERO, hcellmath);
}

void DustFluid::clearDerivedFields(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln)
{
  math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);
  DUST_FLUID_APPLY_TO_DERIVED_FIELDS_ARGS(EXTRA_ARRAY_ZERO, hcellmath);
}

  
void DustFluid::addFieldsToList(std::vector<idx_t> &list)
{
  DUST_FLUID_APPLY_TO_FIELDS_ARGS(ADD_TO_LIST, list);
}

/**
 * @brief  some initialization for each step, currently does nothing
 * 
 */
void DustFluid::stepInit(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{

}

/**
 * @brief  RK evolve physical boundary for particular patch
 * 
 */
void DustFluid::RKEvolvePatchBD(
  const std::shared_ptr<hier::Patch> & patch,
  real_t dt)
{
  std::shared_ptr<hier::PatchGeometry> geom (patch->getPatchGeometry());

  idx_t codim = 1;
  const std::vector<hier::BoundaryBox> & codim1_boxes =
    geom->getCodimensionBoundaries(codim);

  const idx_t n_codim1_boxes = static_cast<idx_t>(codim1_boxes.size());

  // if it has no codimention 1 boundary, it has no other type of boundaries
  if(n_codim1_boxes == 0) return;

  else
    TBOX_ERROR("Current version cannot evolve boudary");
  
}

/**
 * @brief RK evolve patch interior
 */
void DustFluid::RKEvolvePatch(
  const std::shared_ptr<hier::Patch> & patch, BSSN *bssn, real_t dt)
{
  
  // might not need this function
}

void DustFluid::RKEvolvePt(
  idx_t i, idx_t j, idx_t k, BSSNData &bd, DustFluidData & dd, const real_t dx[], real_t dt)
{
  getDustFluidData(i, j, k, &bd, &dd, dx);
  DUST_FLUID_RK_EVOLVE_PT;
}

void DustFluid::RKEvolvePtBd(
  idx_t i, idx_t j, idx_t k, BSSNData &bd, DustFluidData &dd,
  const real_t dx[], real_t dt, int l_idx, int codim)
{
  getDustFluidDataBd(i, j, k, &bd, &dd, dx);
  DUST_FLUID_RK_EVOLVE_BD;
}


void DustFluid::prepareForK1(
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
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(DF_D_a_idx)->getTime()) < EPS)
    {
#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            DUST_FLUID_APPLY_TO_FIELDS(RK4_INIT_R_K1);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(DF_D_p_idx)->getTime())
       - (patch->getPatchData(DF_D_a_idx)->getTime()
          - patch->getPatchData(DF_D_p_idx)->getTime())) < EPS)
    {
          #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            DUST_FLUID_APPLY_TO_FIELDS(RK4_INIT_L_K1);
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
void DustFluid::prepareForK2(
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
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(DF_D_a_idx)->getTime()) < EPS)
    {
          #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            DUST_FLUID_APPLY_TO_FIELDS(RK4_INIT_R_K2);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(DF_D_p_idx)->getTime())
       - (patch->getPatchData(DF_D_a_idx)->getTime()
          - patch->getPatchData(DF_D_p_idx)->getTime())) < EPS)
    {
                #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            DUST_FLUID_APPLY_TO_FIELDS(RK4_INIT_L_K2);
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
void DustFluid::prepareForK3(
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
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(DF_D_a_idx)->getTime()) < EPS)
    {
          #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            DUST_FLUID_APPLY_TO_FIELDS(RK4_INIT_R_K3);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(DF_D_p_idx)->getTime())
       - (patch->getPatchData(DF_D_a_idx)->getTime()
          - patch->getPatchData(DF_D_p_idx)->getTime())) < EPS)
    {
          #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            DUST_FLUID_APPLY_TO_FIELDS(RK4_INIT_L_K3);
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
void DustFluid::prepareForK4(
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
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(DF_D_a_idx)->getTime()) < EPS)
    {
          #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            DUST_FLUID_APPLY_TO_FIELDS(RK4_INIT_R_K4);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(DF_D_p_idx)->getTime())
       - (patch->getPatchData(DF_D_a_idx)->getTime()
          - patch->getPatchData(DF_D_p_idx)->getTime())) < EPS)
    {
          #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            DUST_FLUID_APPLY_TO_FIELDS(RK4_INIT_L_K4);
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
void DustFluid::registerRKRefinerActive(
  xfer::RefineAlgorithm& refiner,
  std::shared_ptr<hier::RefineOperator> &space_refine_op)
{
  DUST_FLUID_APPLY_TO_FIELDS_ARGS(REGISTER_SPACE_REFINE_A, refiner, space_refine_op);
}

  
/**
 * @brief register "scratch" component of the fields to refiner
 * 
 * @param refiner
 * @param refiner operator
 */
void DustFluid::registerRKRefiner(
  xfer::RefineAlgorithm& refiner,
  std::shared_ptr<hier::RefineOperator> &space_refine_op)
{
  DUST_FLUID_APPLY_TO_FIELDS_ARGS(REGISTER_SPACE_REFINE_S, refiner, space_refine_op);
}

/**
 * @brief register "active" component of the fields to coarsener 
 * 
 * @param coarsener
 * @param coarse operator
 */
void DustFluid::registerCoarsenActive(
  xfer::CoarsenAlgorithm& coarsener,
  std::shared_ptr<hier::CoarsenOperator>& coarsen_op)
{
  DUST_FLUID_APPLY_TO_FIELDS_ARGS(REGISTER_COARSEN_A, coarsener, coarsen_op);
}


void DustFluid::copyAToP(
  math::HierarchyCellDataOpsReal<real_t> & hcellmath)
{
  DUST_FLUID_APPLY_TO_FIELDS(COPY_A_TO_P);
}



/**
 * @brief initilizing all pointers for every components's patchdata
 * 
 */
void DustFluid::initPData(
  const std::shared_ptr<hier::Patch> & patch)
{
  DUST_FLUID_APPLY_TO_FIELDS(PDATA_ALL_INIT);
  DUST_FLUID_APPLY_TO_DERIVED_FIELDS_ARGS(PDATA_INIT, a);
}

/**
 * @brief initilizing all pointers for every components's array access
 * 
 */
void DustFluid::initMDA(
  const std::shared_ptr<hier::Patch> & patch)
{
  DUST_FLUID_APPLY_TO_FIELDS(MDA_ACCESS_ALL_INIT);
  DUST_FLUID_APPLY_TO_DERIVED_FIELDS_ARGS(MDA_ACCESS_INIT, a);
}


/**
 * @brief  finalize k1 step for RK on both interior and boundary
 * 
 */
void DustFluid::K1FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);
  
  const hier::Box& box = DF_D_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

      #pragma omp parallel for collapse(2)        
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        DUST_FLUID_APPLY_TO_FIELDS(RK4_FINALIZE_FIELD_1);

      }
    }
  }
}

void DustFluid::K2FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  const hier::Box& box = DF_D_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

    #pragma omp parallel for collapse(2)          
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        DUST_FLUID_APPLY_TO_FIELDS(RK4_FINALIZE_FIELD_2);
      }
    }
  }
}


void DustFluid::K3FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  const hier::Box& box = DF_D_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

      #pragma omp parallel for collapse(2)        
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        DUST_FLUID_APPLY_TO_FIELDS(RK4_FINALIZE_FIELD_3);
      }
    }
  }
}


void DustFluid::K4FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  //Only evolv inner grids for K4
  const hier::Box& box = DF_D_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

    #pragma omp parallel for collapse(2)          
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        DUST_FLUID_APPLY_TO_FIELDS(RK4_FINALIZE_FIELD_4);
      }
    }
  }
}

void DustFluid::printWConstraint(
  BSSN *bssn,   const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t weight_idx)
{
  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  const double * domain_lower = &grid_geometry.getXLower()[0];
  const double * domain_upper = &grid_geometry.getXUpper()[0];

  double L[3];
  for(int i = 0 ; i < DIM; i++)
    L[i] = domain_upper[i] - domain_lower[i];


  real_t max_WV = -1, L2_WV = 0;
  for(idx_t ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
  {
    const std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
    {
      const std::shared_ptr<hier::Patch> & patch = *pit;
      initPData(patch);
      initMDA(patch);

      bssn->initPData(patch);
      bssn->initMDA(patch);

      const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom( 
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      const real_t * dx = &(patch_geom->getDx())[0];
  
      const hier::Box& box = patch->getBox();

      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];

      std::shared_ptr<pdat::CellData<double> > weight(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_idx)));
      

      arr_t weight_array =
        pdat::ArrayDataAccess::access<DIM, double>(
          weight->getArrayData());

      
#pragma omp parallel for collapse(2) reduction(max:max_WV) reduction(+:L2_WV)
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            if(weight_array(i,j,k) > 0)
            {
              BSSNData bd = {0};
              DustFluidData dd = {0};
              bssn->set_bd_values_for_dust_fluid(i, j, k, &bd, dx);
              getDustFluidData(i, j, k, &bd, &dd, dx);
              double W = dd.E / dd.D;
              double v1 = dd.S1 / dd.E, v2 = dd.S2 / dd.E, v3 = dd.S3 / dd.E;
              double Wp = 1.0 / sqrt(1- pw2(bd.chi) * bd.gammai11 * v1 * v1 + pw2(bd.chi) * bd.gammai22 * v2 * v2 + pw2(bd.chi) * bd.gammai33 * v3 * v3
                                         + 2.0 * pw2(bd.chi) * bd.gammai12 * v1 * v2 + 2.0 * pw2(bd.chi) * bd.gammai13 * v1 * v3 + 2.0 * pw2(bd.chi) * bd.gammai23 * v2 * v3);
              max_WV = std::max(max_WV, fabs(Wp));
              L2_WV += pw2(Wp-W)  * weight_array(i,j,k) / (L[0] * L[1] * L[2]);;
            }
          }
        }
      }
    }

  }

  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  mpi.Barrier();
  if (mpi.getSize() > 1) {
    mpi.AllReduce(&max_WV, 1, MPI_MAX);
    mpi.AllReduce(&L2_WV, 1, MPI_SUM);

  }
  L2_WV = sqrt(L2_WV);
  tbox::pout<<"Maximum W is "<<max_WV<<", L2 is "<<L2_WV<<"\n";

}

void DustFluid::addDerivedFields(
  BSSN *bssn,   const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  for(idx_t ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
  {
    const std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));
    addDerivedFields(bssn, level);
  }
}

  
void DustFluid::addDerivedFields(
  BSSN *bssn, const std::shared_ptr<hier::PatchLevel> & level)
{
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
    addDerivedFields(bssn, patch);
  }
}


void DustFluid::addDerivedFields(
  BSSN *bssn, const std::shared_ptr<hier::Patch> & patch )
{
  initPData(patch);
  initMDA(patch);

  bssn->initPData(patch);
  bssn->initMDA(patch);

  const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom( 
    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  const real_t * dx = &(patch_geom->getDx())[0];
  
  const hier::Box& ghost_box = DF_D_a_pdata->getGhostBox();

  hier::Box box = ghost_box;
  
  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

#pragma omp parallel for collapse(2)          
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        BSSNData bd = {0};
        DustFluidData dd = {0};
        bssn->set_bd_values_for_dust_fluid(i, j, k, &bd, dx, false);
        getDustFluidData(i, j, k, &bd, &dd, dx);

        // F01_a(i, j, k) = (bd.alpha * dd.vi1 * dd.D - bd.beta1 * dd.D) / pw3(bd.chi);
        // F02_a(i, j, k) = (bd.alpha * dd.vi2 * dd.D - bd.beta2 * dd.D) / pw3(bd.chi);
        // F03_a(i, j, k) = (bd.alpha * dd.vi3 * dd.D - bd.beta3 * dd.D) / pw3(bd.chi);

        // F11_a(i, j, k) = (bd.alpha *
        //   (dd.m11 * dd.S11 + dd.m12 * dd.S12 + dd.m13 * dd.S13)
        //                   - bd.beta1 * dd.S1) / pw3(bd.chi);
        // F12_a(i, j, k) = (bd.alpha *
        //   (dd.m11 * dd.S12 + dd.m12 * dd.S22 + dd.m13 * dd.S23)
        //                   - bd.beta2 * dd.S1) / pw3(bd.chi);
        // F13_a(i, j, k) = (bd.alpha *
        //   (dd.m11 * dd.S13 + dd.m12 * dd.S23 + dd.m13 * dd.S33)
        //                   - bd.beta3 * dd.S1) / pw3(bd.chi);

        // F21_a(i, j, k) = (bd.alpha *
        //   (dd.m21 * dd.S11 + dd.m22 * dd.S12 + dd.m23 * dd.S13)
        //                   - bd.beta1 * dd.S2) / pw3(bd.chi);
        // F22_a(i, j, k) = (bd.alpha *
        //   (dd.m21 * dd.S12 + dd.m22 * dd.S22 + dd.m23 * dd.S23)
        //                   - bd.beta2 * dd.S2) / pw3(bd.chi);
        // F23_a(i, j, k) = (bd.alpha *
        //   (dd.m31 * dd.S13 + dd.m32 * dd.S23 + dd.m33 * dd.S33)
        //                   - bd.beta3 * dd.S2) / pw3(bd.chi);

        // F31_a(i, j, k) = (bd.alpha *
        //   (dd.m31 * dd.S11 + dd.m32 * dd.S12 + dd.m33 * dd.S13)
        //                   - bd.beta1 * dd.S3) / pw3(bd.chi);
        // F32_a(i, j, k) = (bd.alpha *
        //   (dd.m31 * dd.S12 + dd.m32 * dd.S22 + dd.m33 * dd.S23)
        //                   - bd.beta2 * dd.S3) / pw3(bd.chi);
        // F33_a(i, j, k) = (bd.alpha *
        //   (dd.m31 * dd.S13 + dd.m32 * dd.S23 + dd.m33 * dd.S33)
        //                   - bd.beta3 * dd.S3) / pw3(bd.chi);

        // F41_a(i, j, k) = (bd.alpha * dd.S1 - bd.beta1 * dd.E
        // ) / pw3(bd.chi);
        // F42_a(i, j, k) = (bd.alpha * dd.S2 - bd.beta2 * dd.E
        // ) / pw3(bd.chi);
        // F43_a(i, j, k) = (bd.alpha * dd.S3 - bd.beta3 * dd.E
        // ) / pw3(bd.chi);

        // F01_a(i, j, k) = (bd.alpha * dd.vi1 * dd.D ) ;
        // F02_a(i, j, k) = (bd.alpha * dd.vi2 * dd.D ) ;
        // F03_a(i, j, k) = (bd.alpha * dd.vi3 * dd.D ) ;


        
        // F11_a(i, j, k) = (bd.alpha *
        //                   (dd.m11 * dd.S11 + dd.m12 * dd.S12 + dd.m13 * dd.S13)
        //                   ) ;
        // F12_a(i, j, k) = (bd.alpha *
        //                   (dd.m11 * dd.S12 + dd.m12 * dd.S22 + dd.m13 * dd.S23)
        //                   ) ;
        // F13_a(i, j, k) = (bd.alpha *
        //                   (dd.m11 * dd.S13 + dd.m12 * dd.S23 + dd.m13 * dd.S33)
        //                   ) ;

        // F21_a(i, j, k) = (bd.alpha *
        //                   (dd.m21 * dd.S11 + dd.m22 * dd.S12 + dd.m23 * dd.S13)
        //                   ) ;
        // F22_a(i, j, k) = (bd.alpha *
        //                   (dd.m21 * dd.S12 + dd.m22 * dd.S22 + dd.m23 * dd.S23)
        //                   ) ;
        // F23_a(i, j, k) = (bd.alpha *
        //                   (dd.m21 * dd.S13 + dd.m22 * dd.S23 + dd.m23 * dd.S33)
        //                   ) ;

        // F31_a(i, j, k) = (bd.alpha *
        //                   (dd.m31 * dd.S11 + dd.m32 * dd.S12 + dd.m33 * dd.S13)
        //                   ) ;
        // F32_a(i, j, k) = (bd.alpha *
        //                   (dd.m31 * dd.S12 + dd.m32 * dd.S22 + dd.m33 * dd.S23)
        //                   ) ;
        // F33_a(i, j, k) = (bd.alpha *
        //                   (dd.m31 * dd.S13 + dd.m32 * dd.S23 + dd.m33 * dd.S33)
        //                   ) ;

        F01_a(i, j, k) = (dd.vi1) ;
        F02_a(i, j, k) = (dd.vi2) ;
        F03_a(i, j, k) = (dd.vi3) ;

        
        F11_a(i, j, k) = (bd.alpha * dd.Si1) / dd.E;
        F12_a(i, j, k) = (bd.alpha * dd.Si2) / dd.E;
        F13_a(i, j, k) = (bd.alpha * dd.Si3) / dd.E;

        F21_a(i, j, k) = (bd.alpha * dd.Si1) / dd.E;
        F22_a(i, j, k) = (bd.alpha * dd.Si2) / dd.E;
        F23_a(i, j, k) = (bd.alpha * dd.Si3) / dd.E;

        F31_a(i, j, k) = (bd.alpha * dd.Si1) / dd.E;
        F32_a(i, j, k) = (bd.alpha * dd.Si2) / dd.E;
        F33_a(i, j, k) = (bd.alpha * dd.Si3) / dd.E;

        
        F41_a(i, j, k) = (bd.alpha * dd.Si1 
        ) ;
        F42_a(i, j, k) = (bd.alpha * dd.Si2 
        ) ;
        F43_a(i, j, k) = (bd.alpha * dd.Si3 
        ) ;


        DF_E_a(i, j, k) = dd.E;
      }
    }
  }
  
}

real_t DustFluid::ev_DF_D(BSSNData *bd, DustFluidData *dd, const real_t dx[])
{
#if USE_DUST_FLUID
  if(fabs(dd->D) > 1e30) return 0;
  return
    + 3.0 / bd->chi * (bd->dchidt * dd->D)
    + 3.0 * (bd->d1chi * ( dd->D * bd->alpha * F01_a(bd->i, bd->j, bd->k) - bd->beta1 * dd->D
             )
             + bd->d2chi * (dd->D * bd->alpha * F02_a(bd->i, bd->j, bd->k) - bd->beta2 * dd->D
             )
             + bd->d3chi * (dd->D * bd->alpha * F03_a(bd->i, bd->j, bd->k) - bd->beta3 * dd->D
             )
    ) / bd->chi
    /************** direct scheme *********************************/
    // - derivative(bd->i, bd->j, bd->k, 1, F01_a, dx)
    // - derivative(bd->i, bd->j, bd->k, 2, F02_a, dx)
    // - derivative(bd->i, bd->j, bd->k, 3, F03_a, dx)
    /************************* naive upwind scheme ****************/
    // - upwind_derivative(bd->i, bd->j, bd->k, 1, F01_a, dx, SIGN(bd->x)) * SIGN(bd->x)
    // - upwind_derivative(bd->i, bd->j, bd->k, 2, F02_a, dx, SIGN(bd->y)) * SIGN(bd->y)
    // - upwind_derivative(bd->i, bd->j, bd->k, 3, F03_a, dx, SIGN(bd->z)) * SIGN(bd->z)
    /************************* better upwind scheme? ***************/
    - dd->D * (bd->d1a * dd->vi1 + bd->d2a * dd->vi2 + bd->d3a * dd->vi3)
    - dd->D * bd->alpha * (
      + derivative(bd->i, bd->j, bd->k, 1, F01_a, dx)
      + derivative(bd->i, bd->j, bd->k, 2, F02_a, dx)
      + derivative(bd->i, bd->j, bd->k, 3, F03_a, dx)
    )
    + bd->alpha * (
      + pure_upwind_derivative(bd->i, bd->j, bd->k, 1, DF_D_a, dx, -dd->vi1)
      + pure_upwind_derivative(bd->i, bd->j, bd->k, 2, DF_D_a, dx, -dd->vi2)
      + pure_upwind_derivative(bd->i, bd->j, bd->k, 3, DF_D_a, dx, -dd->vi3)
    )
    + dd->D * (bd->d1beta1 + bd->d2beta2 + bd->d3beta3)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 1, DF_D_a, dx, bd->beta1)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 2, DF_D_a, dx, bd->beta2)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 3, DF_D_a, dx, bd->beta3);
    // + bd->beta1 * derivative(bd->i, bd->j, bd->k, 1, DF_D_a, dx)
    // + bd->beta2 * derivative(bd->i, bd->j, bd->k, 2, DF_D_a, dx)
    // + bd->beta3 * derivative(bd->i, bd->j, bd->k, 3, DF_D_a, dx);
#endif
}

real_t DustFluid::ev_DF_S1(BSSNData *bd, DustFluidData *dd, const real_t dx[])
{
  // if( bd->i <= 66 && bd->i >= 61 && bd->j <= 66 && bd->j >= 61 && bd->k <= 66 && bd->k >= 61)
  //   return 0;
#if USE_DUST_FLUID
  if(fabs(dd->S1) > 1e30) return 0;
  return 
    (0.5 * bd->alpha *
     (dd->S11 * dd->d1m11 + dd->S22 * dd->d1m22 + dd->S33 * dd->d1m33
      + 2.0 * dd->S12 * dd->d1m12 + 2.0 * dd->S13 * dd->d1m13 + 2.0 * dd->S23 * dd->d1m23)
     + (dd->S1 * bd->d1beta1 + dd->S2 * bd->d1beta2 + dd->S3 * bd->d1beta3)
     - dd->E * bd->d1a
    )
    + 3.0 / bd->chi * (bd->dchidt * dd->S1)
    + 3.0 * (bd->d1chi * ( F11_a(bd->i, bd->j, bd->k) * dd->S1- bd->beta1 * dd->S1
             )
             + bd->d2chi * (F12_a(bd->i, bd->j, bd->k) * dd->S1- bd->beta2 * dd->S1
             )
             + bd->d3chi * (F13_a(bd->i, bd->j, bd->k) * dd->S1- bd->beta3 * dd->S1
             )
    ) / bd->chi
    /************** direct scheme *********************************/
    // - derivative(bd->i, bd->j, bd->k, 1, F11_a, dx)
    // - derivative(bd->i, bd->j, bd->k, 2, F12_a, dx)
    // - derivative(bd->i, bd->j, bd->k, 3, F13_a, dx)
    /************************* naive upwind scheme ****************/
    // - upwind_derivative(bd->i, bd->j, bd->k, 1, F11_a, dx, SIGN(bd->x)) * SIGN(bd->x)
    // - upwind_derivative(bd->i, bd->j, bd->k, 2, F12_a, dx, SIGN(bd->y)) * SIGN(bd->y)
    // - upwind_derivative(bd->i, bd->j, bd->k, 3, F13_a, dx, SIGN(bd->z)) * SIGN(bd->z)
    /************************* better upwind scheme? ***************/
    - dd->S1 * (
      + derivative(bd->i, bd->j, bd->k, 1, F11_a, dx)
      + derivative(bd->i, bd->j, bd->k, 2, F12_a, dx)
      + derivative(bd->i, bd->j, bd->k, 3, F13_a, dx))
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 1, DF_S1_a, dx, -F11_a(bd->i, bd->j,bd->k))
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 2, DF_S1_a, dx, -F12_a(bd->i, bd->j,bd->k))
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 3, DF_S1_a, dx, -F13_a(bd->i, bd->j,bd->k))
    + dd->S1 * (bd->d1beta1 + bd->d2beta2 + bd->d3beta3)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 1, DF_S1_a, dx, bd->beta1)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 2, DF_S1_a, dx, bd->beta2)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 3, DF_S1_a, dx, bd->beta3);
    // + bd->beta1 * derivative(bd->i, bd->j, bd->k, 1, DF_S1_a, dx)
    // + bd->beta2 * derivative(bd->i, bd->j, bd->k, 2, DF_S1_a, dx)
    // + bd->beta3 * derivative(bd->i, bd->j, bd->k, 3, DF_S1_a, dx);

#endif
}

real_t DustFluid::ev_DF_S2(BSSNData *bd, DustFluidData *dd, const real_t dx[])
{
  // if( bd->i <= 66 && bd->i >= 61 && bd->j <= 66 && bd->j >= 61 && bd->k <= 66 && bd->k >= 61)
  //   return 0;
#if USE_DUST_FLUID
  if(fabs(dd->S2) > 1e30) return 0;
  return 
    (0.5 * bd->alpha *
     (dd->S11 * dd->d2m11 + dd->S22 * dd->d2m22 + dd->S33 * dd->d2m33
      + 2.0 * dd->S12 * dd->d2m12 + 2.0 * dd->S13 * dd->d2m13 + 2.0 * dd->S23 * dd->d2m23)
     + (dd->S1 * bd->d2beta1 + dd->S2 * bd->d2beta2 + dd->S3 * bd->d2beta3)
     - dd->E * bd->d2a
    )
    + 3.0 / bd->chi * (bd->dchidt * dd->S2)
    + 3.0 * (bd->d1chi * ( F21_a(bd->i, bd->j, bd->k) * dd->S2 - bd->beta1 * dd->S2
             )
             + bd->d2chi * (F22_a(bd->i, bd->j, bd->k) * dd->S2 - bd->beta2 * dd->S2
             )
             + bd->d3chi * (F23_a(bd->i, bd->j, bd->k) * dd->S2 - bd->beta3 * dd->S2
             )
    ) / bd->chi
    // - derivative(bd->i, bd->j, bd->k, 1, F21_a, dx)
    // - derivative(bd->i, bd->j, bd->k, 2, F22_a, dx)
    // - derivative(bd->i, bd->j, bd->k, 3, F23_a, dx)
    // - upwind_derivative(bd->i, bd->j, bd->k, 1, F21_a, dx, SIGN(bd->x)) * SIGN(bd->x)
    // - upwind_derivative(bd->i, bd->j, bd->k, 2, F22_a, dx, SIGN(bd->y)) * SIGN(bd->y)
    // - upwind_derivative(bd->i, bd->j, bd->k, 3, F23_a, dx, SIGN(bd->z)) * SIGN(bd->z)
    - dd->S2 * (
      + derivative(bd->i, bd->j, bd->k, 1, F21_a, dx)
      + derivative(bd->i, bd->j, bd->k, 2, F22_a, dx)
      + derivative(bd->i, bd->j, bd->k, 3, F23_a, dx))
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 1, DF_S2_a, dx, -F21_a(bd->i, bd->j,bd->k))
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 2, DF_S2_a, dx, -F22_a(bd->i, bd->j,bd->k))
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 3, DF_S2_a, dx, -F23_a(bd->i, bd->j,bd->k))

    + dd->S2 * (bd->d1beta1 + bd->d2beta2 + bd->d3beta3)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 1, DF_S2_a, dx, bd->beta1)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 2, DF_S2_a, dx, bd->beta2)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 3, DF_S2_a, dx, bd->beta3);
    // + bd->beta1 * derivative(bd->i, bd->j, bd->k, 1, DF_S2_a, dx)
    // + bd->beta2 * derivative(bd->i, bd->j, bd->k, 2, DF_S2_a, dx)
    // + bd->beta3 * derivative(bd->i, bd->j, bd->k, 3, DF_S2_a, dx);
#endif

}

real_t DustFluid::ev_DF_S3(BSSNData *bd, DustFluidData *dd, const real_t dx[])
{
  // if( bd->i <= 66 && bd->i >= 61 && bd->j <= 66 && bd->j >= 61 && bd->k <= 66 && bd->k >= 61)
  //   return 0;
#if USE_DUST_FLUID
  if(fabs(dd->S3) > 1e30) return 0;
  return 
    (0.5 * bd->alpha *
     (dd->S11 * dd->d3m11 + dd->S22 * dd->d3m22 + dd->S33 * dd->d3m33
      + 2.0 * dd->S12 * dd->d3m12 + 2.0 * dd->S13 * dd->d3m13 + 2.0 * dd->S23 * dd->d3m23)
     + (dd->S1 * bd->d3beta1 + dd->S2 * bd->d3beta2 + dd->S3 * bd->d3beta3)
     - dd->E * bd->d3a
    )
    + 3.0 / bd->chi * (bd->dchidt * dd->S3)
    + 3.0 * (bd->d1chi * ( F31_a(bd->i, bd->j, bd->k) * dd->S3 - bd->beta1 * dd->S3
             )
             + bd->d2chi * (F32_a(bd->i, bd->j, bd->k) * dd->S3 - bd->beta2 * dd->S3
             )
             + bd->d3chi * (F33_a(bd->i, bd->j, bd->k) * dd->S3 - bd->beta3 * dd->S3
             )
    ) / bd->chi
    // - derivative(bd->i, bd->j, bd->k, 1, F31_a, dx)
    // - derivative(bd->i, bd->j, bd->k, 2, F32_a, dx)
    // - derivative(bd->i, bd->j, bd->k, 3, F33_a, dx)
    // - upwind_derivative(bd->i, bd->j, bd->k, 1, F31_a, dx, SIGN(bd->x)) * SIGN(bd->x)
    // - upwind_derivative(bd->i, bd->j, bd->k, 2, F32_a, dx, SIGN(bd->y)) * SIGN(bd->y)
    // - upwind_derivative(bd->i, bd->j, bd->k, 3, F33_a, dx, SIGN(bd->z)) * SIGN(bd->z)
    - dd->S3 * (
      + derivative(bd->i, bd->j, bd->k, 1, F31_a, dx)
      + derivative(bd->i, bd->j, bd->k, 2, F32_a, dx)
      + derivative(bd->i, bd->j, bd->k, 3, F33_a, dx))
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 1, DF_S3_a, dx, -F31_a(bd->i, bd->j,bd->k))
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 2, DF_S3_a, dx, -F32_a(bd->i, bd->j,bd->k))
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 3, DF_S3_a, dx, -F33_a(bd->i, bd->j,bd->k))

    + dd->S3* (bd->d1beta1 + bd->d2beta2 + bd->d3beta3)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 1, DF_S3_a, dx, bd->beta1)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 2, DF_S3_a, dx, bd->beta2)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 3, DF_S3_a, dx, bd->beta3);
    // + bd->beta1 * derivative(bd->i, bd->j, bd->k, 1, DF_S3_a, dx)
    // + bd->beta2 * derivative(bd->i, bd->j, bd->k, 2, DF_S3_a, dx)
    // + bd->beta3 * derivative(bd->i, bd->j, bd->k, 3, DF_S3_a, dx);

#endif
}

real_t DustFluid::ev_DF_E(BSSNData *bd, DustFluidData *dd, const real_t dx[])
{
#if USE_DUST_FLUID
  return 0;
  return
    (bd->alpha *
     (dd->S11 * dd->K11  + dd->S22 * dd->K22 + dd->S33 * dd->K33
      + 2.0 * dd->S12 * dd->K12 + 2.0 * dd->S13 * dd->K13 + 2.0 * dd->S23 * dd->K23)
     - (dd->Si1 * bd->d1a + dd->Si2 * bd->d2a + dd->Si3 * bd->d3a)
    )/ pw3(bd->chi); 
  - derivative(bd->i, bd->j, bd->k, 1, F41_a, dx)
    - derivative(bd->i, bd->j, bd->k, 2, F42_a, dx)
    - derivative(bd->i, bd->j, bd->k, 3, F43_a, dx)
    + DF_E_a(bd->i, bd->j, bd->k) * (bd->d1beta1 + bd->d2beta2 + bd->d3beta3)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 1, DF_E_a, dx, bd->beta1)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 2, DF_E_a, dx, bd->beta2)
    + pure_upwind_derivative(bd->i, bd->j, bd->k, 3, DF_E_a, dx, bd->beta3);
    // + bd->beta1 * derivative(bd->i, bd->j, bd->k, 1, DF_E_a, dx)
    // + bd->beta2 * derivative(bd->i, bd->j, bd->k, 2, DF_E_a, dx)
    // + bd->beta3 * derivative(bd->i, bd->j, bd->k, 3, DF_E_a, dx);
#endif
}

real_t DustFluid::ev_DF_D_bd(BSSNData *bd, DustFluidData *dd, const real_t dx[], int l_idx, int codim)
{
  
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, DF_D_a, dx, l_idx, codim) * bd->x
                        + bd_derivative(bd->i, bd->j, bd->k, 2, DF_D_a, dx, l_idx, codim) * bd->y
                        + bd_derivative(bd->i, bd->j, bd->k, 3, DF_D_a, dx, l_idx, codim) * bd->z 
                        + dd->D           );

}

real_t DustFluid::ev_DF_S1_bd(BSSNData *bd, DustFluidData *dd, const real_t dx[], int l_idx, int codim)
{
    return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, DF_S1_a, dx, l_idx, codim) * bd->x
                        + bd_derivative(bd->i, bd->j, bd->k, 2, DF_S1_a, dx, l_idx, codim) * bd->y
                        + bd_derivative(bd->i, bd->j, bd->k, 3, DF_S1_a, dx, l_idx, codim) * bd->z 
                        + dd->S1           );
}

real_t DustFluid::ev_DF_S2_bd(BSSNData *bd, DustFluidData *dd, const real_t dx[], int l_idx, int codim)
{
    return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, DF_S2_a, dx, l_idx, codim) * bd->x
                        + bd_derivative(bd->i, bd->j, bd->k, 2, DF_S2_a, dx, l_idx, codim) * bd->y
                        + bd_derivative(bd->i, bd->j, bd->k, 3, DF_S2_a, dx, l_idx, codim) * bd->z 
                        + dd->S2           );
}

real_t DustFluid::ev_DF_S3_bd(BSSNData *bd, DustFluidData *dd, const real_t dx[], int l_idx, int codim)
{
    return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, DF_S3_a, dx, l_idx, codim) * bd->x
                        + bd_derivative(bd->i, bd->j, bd->k, 2, DF_S3_a, dx, l_idx, codim) * bd->y
                        + bd_derivative(bd->i, bd->j, bd->k, 3, DF_S3_a, dx, l_idx, codim) * bd->z 
                        + dd->S3           );
}

real_t DustFluid::ev_DF_E_bd(BSSNData *bd, DustFluidData *dd, const real_t dx[], int l_idx, int codim)
{
    return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, DF_E_a, dx, l_idx, codim) * bd->x
                        + bd_derivative(bd->i, bd->j, bd->k, 2, DF_E_a, dx, l_idx, codim) * bd->y
                        + bd_derivative(bd->i, bd->j, bd->k, 3, DF_E_a, dx, l_idx, codim) * bd->z 
                        + dd->E           );

}

void DustFluid::getDustFluidData(
  idx_t i, idx_t j, idx_t k, BSSNData *bd, DustFluidData * dd, const real_t dx[])
{
  dd->D = DF_D_a(i, j, k) ;
  dd->S1 = DF_S1_a(i, j, k) ;
  dd->S2 = DF_S2_a(i, j, k) ;
  dd->S3 = DF_S3_a(i, j, k) ;
  //  dd->E = DF_E_a(i, j, k) * pw3(bd->chi);

  double Wp =  sqrt(1.0 + pw2(bd->chi) * (bd->gammai11 * dd->S1 * dd->S1
                    + bd->gammai22 * dd->S2 * dd->S2
                    + bd->gammai33 * dd->S3 * dd->S3
                    + 2.0 * bd->gammai12 * dd->S1 * dd->S2
                    + 2.0 * bd->gammai13 * dd->S1 * dd->S3
                    + 2.0 * bd->gammai23 * dd->S2 * dd->S3) / pw2(dd->D));

  dd->E = dd->D * Wp;


  dd->v1 = dd->S1 / dd->E;
  dd->v2 = dd->S2 / dd->E;
  dd->v3 = dd->S3 / dd->E;


  
  dd->m11 = bd->gamma11 / pw2(bd->chi);
  dd->m22 = bd->gamma22 / pw2(bd->chi);
  dd->m22 = bd->gamma33 / pw2(bd->chi);
  dd->m12 = bd->gamma12 / pw2(bd->chi);
  dd->m13 = bd->gamma13 / pw2(bd->chi);
  dd->m23 = bd->gamma23 / pw2(bd->chi);

  dd->K11 = (bd->A11 + bd->gamma11 * bd->K / 3.0) / pw2(bd->chi);
  dd->K22 = (bd->A22 + bd->gamma22 * bd->K / 3.0) / pw2(bd->chi);
  dd->K33 = (bd->A33 + bd->gamma33 * bd->K / 3.0) / pw2(bd->chi);
  dd->K12 = (bd->A12 + bd->gamma12 * bd->K / 3.0) / pw2(bd->chi);
  dd->K13 = (bd->A13 + bd->gamma13 * bd->K / 3.0) / pw2(bd->chi);
  dd->K23 = (bd->A23 + bd->gamma23 * bd->K / 3.0) / pw2(bd->chi);
  
  dd->vi1 = pw2(bd->chi) * bd->gammai11 * dd->v1
    + pw2(bd->chi) * bd->gammai12 * dd->v2
    + pw2(bd->chi) * bd->gammai13 * dd->v3;
  dd->vi2 = pw2(bd->chi) * bd->gammai21 * dd->v1
    + pw2(bd->chi) * bd->gammai22 * dd->v2
    + pw2(bd->chi) * bd->gammai23 * dd->v3;
  dd->vi3 = pw2(bd->chi) * bd->gammai31 * dd->v1
    + pw2(bd->chi) * bd->gammai32 * dd->v2
    + pw2(bd->chi) * bd->gammai33 * dd->v3;

  dd->Si1 = dd->vi1 * dd->E;
  dd->Si2 = dd->vi2 * dd->E;
  dd->Si3 = dd->vi3 * dd->E;

  dd->S11 = dd->E * dd->vi1 * dd->vi1;
  dd->S22 = dd->E * dd->vi2 * dd->vi2;
  dd->S33 = dd->E * dd->vi3 * dd->vi3;
  dd->S12 = dd->E * dd->vi1 * dd->vi2;
  dd->S13 = dd->E * dd->vi1 * dd->vi3;
  dd->S23 = dd->E * dd->vi2 * dd->vi3;


  
  dd->d1m11 = -2.0 * bd->d1chi / pw3(bd->chi) * dd->m11 + bd->d1g11 / pw2(bd->chi) ;
  dd->d1m22 = -2.0 * bd->d1chi / pw3(bd->chi) * dd->m22 + bd->d1g22 / pw2(bd->chi) ;
  dd->d1m33 = -2.0 * bd->d1chi / pw3(bd->chi) * dd->m33 + bd->d1g33 / pw2(bd->chi) ; 
  dd->d1m12 = -2.0 * bd->d1chi / pw3(bd->chi) * dd->m12 + bd->d1g12 / pw2(bd->chi) ;
  dd->d1m13 = -2.0 * bd->d1chi / pw3(bd->chi) * dd->m13 + bd->d1g13 / pw2(bd->chi) ;
  dd->d1m23 = -2.0 * bd->d1chi / pw3(bd->chi) * dd->m23 + bd->d1g23 / pw2(bd->chi) ;

  dd->d2m11 = -2.0 * bd->d2chi / pw3(bd->chi) * dd->m11 + bd->d2g11 / pw2(bd->chi) ;
  dd->d2m22 = -2.0 * bd->d2chi / pw3(bd->chi) * dd->m22 + bd->d2g22 / pw2(bd->chi) ;
  dd->d2m33 = -2.0 * bd->d2chi / pw3(bd->chi) * dd->m33 + bd->d2g33 / pw2(bd->chi) ; 
  dd->d2m12 = -2.0 * bd->d2chi / pw3(bd->chi) * dd->m12 + bd->d2g12 / pw2(bd->chi) ;
  dd->d2m13 = -2.0 * bd->d2chi / pw3(bd->chi) * dd->m13 + bd->d2g13 / pw2(bd->chi) ;
  dd->d2m23 = -2.0 * bd->d2chi / pw3(bd->chi) * dd->m23 + bd->d2g23 / pw2(bd->chi) ;

  dd->d3m11 = -2.0 * bd->d3chi / pw3(bd->chi) * dd->m11 + bd->d3g11 / pw2(bd->chi) ;
  dd->d3m22 = -2.0 * bd->d3chi / pw3(bd->chi) * dd->m22 + bd->d3g22 / pw2(bd->chi) ;
  dd->d3m33 = -2.0 * bd->d3chi / pw3(bd->chi) * dd->m33 + bd->d3g33 / pw2(bd->chi) ; 
  dd->d3m12 = -2.0 * bd->d3chi / pw3(bd->chi) * dd->m12 + bd->d3g12 / pw2(bd->chi) ;
  dd->d3m13 = -2.0 * bd->d3chi / pw3(bd->chi) * dd->m13 + bd->d3g13 / pw2(bd->chi) ;
  dd->d3m23 = -2.0 * bd->d3chi / pw3(bd->chi) * dd->m23 + bd->d3g23 / pw2(bd->chi) ;

}
    
void DustFluid::getDustFluidDataBd(
  idx_t i, idx_t j, idx_t k, BSSNData *bd, DustFluidData * dd, const real_t dx[])
{
#if USE_SOMMERFIELD_BOUNDARY
  dd->D = DF_D_a(i, j, k) * pw3(bd->chi);
  dd->S1 = DF_S1_a(i, j, k) * pw3(bd->chi) ;
  dd->S2 = DF_S2_a(i, j, k) * pw3(bd->chi) ;
  dd->S3 = DF_S3_a(i, j, k) * pw3(bd->chi);
  dd->E = DF_E_a(i, j, k) * pw3(bd->chi);
#endif
}


void DustFluid::addBSSNSrc(
  BSSN *bssn,   const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  for(idx_t ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
  {
    const std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));
    addBSSNSrc(bssn, level);
  }
}

  
void DustFluid::addBSSNSrc(
  BSSN *bssn, const std::shared_ptr<hier::PatchLevel> & level)
{
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
    addBSSNSrc(bssn, patch, true);
  }
}


  
void DustFluid::addBSSNSrc(
  BSSN *bssn, const std::shared_ptr<hier::Patch> & patch, bool need_init_arr)
{
  if(is_test_fluid == true)
    return;
  else
    TBOX_ERROR("Non test fluid feature is not ready yet!") ;
  if(need_init_arr)
  {
    bssn->initPData(patch);
    bssn->initMDA(patch);

    initPData(patch);
    initMDA(patch);
  }

  arr_t & DIFFr_a = bssn->DIFFr_a;
  arr_t & DIFFS_a = bssn->DIFFS_a;
  arr_t & BSSN_S1_a = bssn->S1_a;
  arr_t & BSSN_S2_a = bssn->S2_a;
  arr_t & BSSN_S3_a = bssn->S3_a;
  arr_t & STF11_a = bssn->STF11_a;
  arr_t & STF12_a = bssn->STF12_a;
  arr_t & STF13_a = bssn->STF13_a;
  arr_t & STF22_a = bssn->STF22_a;
  arr_t & STF23_a = bssn->STF23_a;
  arr_t & STF33_a = bssn->STF33_a;

  const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom( 
    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

  const real_t * dx = &(patch_geom->getDx())[0];

  const hier::Box& box = patch->getBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

#pragma omp parallel for collapse(2)        
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        BSSNData bd = {0};
        DustFluidData dd = {0};
        // TODO: remove redundant computations here?
        bssn->set_bd_values_for_dust_fluid(i, j, k, &bd, dx);
        // including source data
        getDustFluidData(i, j, k, &bd, &dd, dx);
        
      {
        DIFFr_a(i, j, k) += dd.E;
      }
        
      }
    }
  }
  return;    

  
}


void DustFluid::setLevelTime(
  const std::shared_ptr<hier::PatchLevel> & level,
  double from_t, double to_t)
{
  DUST_FLUID_APPLY_TO_FIELDS_ARGS(SET_LEVEL_TIME, from_t, to_t);
}


} // namespace cosmo
