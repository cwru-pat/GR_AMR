#include "scalar.h"
#include "../../cosmo_includes.h"
#include "../../utils/math.h"


using namespace SAMRAI;

namespace cosmo
{

Scalar::Scalar(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  const tbox::Dimension& dim_in,
  std::shared_ptr<tbox::Database> database_in,
  std::ostream* l_stream_in,
  real_t KO_damping_coefficient_in):
  lstream(l_stream_in),
  cosmo_scalar_db(database_in),
  dim(dim_in),
  KO_damping_coefficient(KO_damping_coefficient_in),
  potentialHandler(new scalarPotentialHandler(cosmo_scalar_db))
{
  SCALAR_APPLY_TO_FIELDS(VAR_INIT);

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
#if USE_BACKUP_FIELDS
  std::shared_ptr<hier::VariableContext> context_b(
    variable_db->getContext("BACKUP"));
#endif

  
  SCALAR_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_scratch, s, STENCIL_ORDER);
  SCALAR_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_previous, p, STENCIL_ORDER);
  SCALAR_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_active, a, STENCIL_ORDER);
  SCALAR_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k1, k1, STENCIL_ORDER);
  SCALAR_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k2, k2, STENCIL_ORDER);
  SCALAR_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k3, k3, STENCIL_ORDER);
  SCALAR_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_k4, k4, STENCIL_ORDER);
#if USE_BACKUP_FIELDS
  SCALAR_APPLY_TO_FIELDS_ARGS(REG_TO_CONTEXT, context_b, b, STENCIL_ORDER);
#endif
    
}

Scalar::~Scalar()
{
  
}

/**
 * @brief  setting length of physical domain and chi lower bound
 * 
 */
void Scalar::init(const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  return;
}
  
void Scalar::alloc(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln)
{
  std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
  SCALAR_APPLY_TO_FIELDS(RK4_ARRAY_ALLOC);
}

void Scalar::clear(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy, idx_t ln)
{
  math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);
  SCALAR_APPLY_TO_FIELDS_ARGS(RK4_ARRAY_ZERO, hcellmath);
}

void Scalar::addFieldsToList(std::vector<idx_t> &list)
{
  SCALAR_APPLY_TO_FIELDS_ARGS(ADD_TO_LIST, list);
}

/**
 * @brief  some initialization for each step, currently does nothing
 * 
 */
void Scalar::stepInit(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
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
void Scalar::RKEvolvePatchBD(
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
void Scalar::RKEvolvePatch(
  const std::shared_ptr<hier::Patch> & patch, real_t dt)
{
  // might not need this function
}

  
void Scalar::RKEvolvePt(
  idx_t i, idx_t j, idx_t k, BSSNData &bd, ScalarData & sd, const real_t dx[], real_t dt)
{
  getScalarData(i, j, k, &bd, &sd, dx);
  SCALAR_RK_EVOLVE_PT;
}

void Scalar::RKEvolvePtBd(
  idx_t i, idx_t j, idx_t k, BSSNData &bd, ScalarData &sd,
  const real_t dx[], real_t dt, int l_idx, int codim)
{
  getScalarDataBd(i, j, k, &bd, &sd, dx);
  SCALAR_RK_EVOLVE_BD;
}

  
void Scalar::prepareForK1(
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
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(phi_a_idx)->getTime()) < EPS)
    {
#pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            SCALAR_APPLY_TO_FIELDS(RK4_INIT_R_K1);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(phi_p_idx)->getTime())
       - (patch->getPatchData(phi_a_idx)->getTime()
          - patch->getPatchData(phi_p_idx)->getTime())) < EPS)
    {
          #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            SCALAR_APPLY_TO_FIELDS(RK4_INIT_L_K1);
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
void Scalar::prepareForK2(
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
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(phi_a_idx)->getTime()) < EPS)
    {
          #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            SCALAR_APPLY_TO_FIELDS(RK4_INIT_R_K2);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(phi_p_idx)->getTime())
       - (patch->getPatchData(phi_a_idx)->getTime()
          - patch->getPatchData(phi_p_idx)->getTime())) < EPS)
    {
                #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            SCALAR_APPLY_TO_FIELDS(RK4_INIT_L_K2);
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
void Scalar::prepareForK3(
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
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(phi_a_idx)->getTime()) < EPS)
    {
          #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            SCALAR_APPLY_TO_FIELDS(RK4_INIT_R_K3);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(phi_p_idx)->getTime())
       - (patch->getPatchData(phi_a_idx)->getTime()
          - patch->getPatchData(phi_p_idx)->getTime())) < EPS)
    {
          #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            SCALAR_APPLY_TO_FIELDS(RK4_INIT_L_K3);
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
void Scalar::prepareForK4(
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
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(phi_a_idx)->getTime()) < EPS)
    {
          #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            SCALAR_APPLY_TO_FIELDS(RK4_INIT_R_K4);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(phi_p_idx)->getTime())
       - (patch->getPatchData(phi_a_idx)->getTime()
          - patch->getPatchData(phi_p_idx)->getTime())) < EPS)
    {
          #pragma omp parallel for collapse(2)        
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            SCALAR_APPLY_TO_FIELDS(RK4_INIT_L_K4);
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
void Scalar::registerRKRefinerActive(
  xfer::RefineAlgorithm& refiner,
  std::shared_ptr<hier::RefineOperator> &space_refine_op)
{
  SCALAR_APPLY_TO_FIELDS_ARGS(REGISTER_SPACE_REFINE_A, refiner, space_refine_op);
}

  
/**
 * @brief register "scratch" component of the fields to refiner
 * 
 * @param refiner
 * @param refiner operator
 */
void Scalar::registerRKRefiner(
  xfer::RefineAlgorithm& refiner,
  std::shared_ptr<hier::RefineOperator> &space_refine_op)
{
  SCALAR_APPLY_TO_FIELDS_ARGS(REGISTER_SPACE_REFINE_S, refiner, space_refine_op);
}

/**
 * @brief register "active" component of the fields to coarsener 
 * 
 * @param coarsener
 * @param coarse operator
 */
void Scalar::registerCoarsenActive(
  xfer::CoarsenAlgorithm& coarsener,
  std::shared_ptr<hier::CoarsenOperator>& coarsen_op)
{
  SCALAR_APPLY_TO_FIELDS_ARGS(REGISTER_COARSEN_A, coarsener, coarsen_op);
}


void Scalar::copyAToP(
  math::HierarchyCellDataOpsReal<real_t> & hcellmath)
{
  SCALAR_APPLY_TO_FIELDS(COPY_A_TO_P);
}

#if USE_BACKUP_FIELDS

void Scalar::copyPToB(
  math::HierarchyCellDataOpsReal<real_t> & hcellmath)
{
  SCALAR_APPLY_TO_FIELDS(COPY_P_TO_B);
}

void Scalar::copyBToP(
  math::HierarchyCellDataOpsReal<real_t> & hcellmath)
{
  SCALAR_APPLY_TO_FIELDS(COPY_B_TO_P);
}

void Scalar::copyBToA(
  math::HierarchyCellDataOpsReal<real_t> & hcellmath)
{
  SCALAR_APPLY_TO_FIELDS(COPY_B_TO_A);
}


#endif




/**
 * @brief initilizing all pointers for every components's patchdata
 * 
 */
void Scalar::initPData(
  const std::shared_ptr<hier::Patch> & patch)
{
  SCALAR_APPLY_TO_FIELDS(PDATA_ALL_INIT);
}

/**
 * @brief initilizing all pointers for every components's array access
 * 
 */
void Scalar::initMDA(
  const std::shared_ptr<hier::Patch> & patch)
{
  SCALAR_APPLY_TO_FIELDS(MDA_ACCESS_ALL_INIT);
}


/**
 * @brief  finalize k1 step for RK on both interior and boundary
 * 
 */
void Scalar::K1FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);
  
  const hier::Box& box = phi_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

      #pragma omp parallel for collapse(2)        
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        SCALAR_APPLY_TO_FIELDS(RK4_FINALIZE_FIELD_1);
      }
    }
  }
}

void Scalar::K2FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  const hier::Box& box = phi_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

    #pragma omp parallel for collapse(2)          
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        SCALAR_APPLY_TO_FIELDS(RK4_FINALIZE_FIELD_2);
      }
    }
  }
}


void Scalar::K3FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  const hier::Box& box = phi_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

      #pragma omp parallel for collapse(2)        
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        SCALAR_APPLY_TO_FIELDS(RK4_FINALIZE_FIELD_3);
      }
    }
  }
}


void Scalar::K4FinalizePatch(
  const std::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  //Only evolv inner grids for K4
  const hier::Box& box = phi_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

    #pragma omp parallel for collapse(2)          
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        SCALAR_APPLY_TO_FIELDS(RK4_FINALIZE_FIELD_4);
      }
    }
  }
}


void Scalar::getScalarData(
  idx_t i, idx_t j, idx_t k, BSSNData *bd, ScalarData * sd, const real_t dx[])
{
  sd->phi = phi_a(i, j, k);
  sd->Pi = Pi_a(i, j, k);

  sd->d1phi = derivative(i, j, k, 1, phi_a, dx);
  sd->d2phi = derivative(i, j, k, 2, phi_a, dx);
  sd->d3phi = derivative(i, j, k, 3, phi_a, dx);
}

void Scalar::getScalarDataBd(
  idx_t i, idx_t j, idx_t k, BSSNData *bd, ScalarData * sd, const real_t dx[])
{
#if USE_SOMMERFIELD_BOUNDARY
  sd->phi = phi_a(i, j, k);
  sd->Pi = Pi_a(i, j, k);
#endif
}


real_t Scalar::ev_phi(BSSNData *bd, ScalarData *sd, const real_t dx[])
{
  return (
    sd->Pi
  )
    - KO_dissipation_Q(bd->i, bd->j, bd->k, phi_a, dx, KO_damping_coefficient);
}

real_t Scalar::ev_Pi(BSSNData *bd, ScalarData *sd, const real_t dx[])
{
  return (
    -3.0 * bd->H * sd->Pi +
    laplacian(bd->i, bd->j, bd->k, phi_a, dx) / pw2(bd->a)
    - potentialHandler->ev_der_potential(bd, sd)
  )
    - KO_dissipation_Q(bd->i, bd->j, bd->k, Pi_a, dx, KO_damping_coefficient);
}


real_t Scalar::ev_phi_bd(BSSNData *bd, ScalarData *sd, const real_t dx[], int l_idx, int codim)
{
  return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, phi_a, dx, l_idx, codim) * bd->x
                        + bd_derivative(bd->i, bd->j, bd->k, 2, phi_a, dx, l_idx, codim) * bd->y
                        + bd_derivative(bd->i, bd->j, bd->k, 3, phi_a, dx, l_idx, codim) * bd->z 
                        + sd->phi           );

}

real_t Scalar::ev_Pi_bd(BSSNData *bd, ScalarData *sd, const real_t dx[], int l_idx, int codim)
{
    return -1.0/bd->norm*(bd_derivative(bd->i, bd->j, bd->k, 1, Pi_a, dx, l_idx, codim) * bd->x
                        + bd_derivative(bd->i, bd->j, bd->k, 2, Pi_a, dx, l_idx, codim) * bd->y
                        + bd_derivative(bd->i, bd->j, bd->k, 3, Pi_a, dx, l_idx, codim) * bd->z 
                        + sd->Pi           );
}
  
void Scalar::addBSSNSrc(
  BSSN * bssn, const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  for(idx_t ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
  {
    const std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));
    addBSSNSrc(bssn, level);
  }
  
}

void Scalar::addBSSNSrc(
  BSSN * bssn, const std::shared_ptr<hier::PatchLevel> & level)
{
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
    addBSSNSrc(bssn, patch, true);
  }  
}

/**
 * @brief set the time of previous components of all BSSN fields as from_t
 *        set the time of active components of all BSSN fields as to_t
 */ 
void Scalar::setLevelTime(
  const std::shared_ptr<hier::PatchLevel> & level,
  double from_t, double to_t)
{
  SCALAR_APPLY_TO_FIELDS_ARGS(SET_LEVEL_TIME, from_t, to_t);
}

  
void Scalar::addBSSNSrc(
  BSSN * bssn, const std::shared_ptr<hier::Patch> & patch, bool need_init_arr)
{
  if(need_init_arr)
  {
    bssn->initPData(patch);
    bssn->initMDA(patch);

    initPData(patch);
    initMDA(patch);
  }
  
  arr_t & DIFFr_a = bssn->DIFFr_a;
  arr_t & DIFFS_a = bssn->DIFFS_a;

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
        ScalarData sd = {0};
        // TODO: remove redundant computations here?
        bssn->set_bd_values(i, j, k, &bd, dx);
        getScalarData(i, j, k, &bd, &sd, dx);


        DIFFr_a(i,j,k) += 0.5*pw2(ev_phi(&bd, &sd, dx))
          + 0.5*(pw2(sd.d1phi)+pw2(sd.d2phi) + pw2(sd.d3phi)) / pw2(bd.a) + potentialHandler->ev_potential(&bd, &sd);

        DIFFS_a(i,j,k) += 0.5*pw2(ev_phi(&bd, &sd, dx))
          - (pw2(sd.d1phi)+pw2(sd.d2phi) + pw2(sd.d3phi)) / pw2(bd.a)/6.0 - potentialHandler->ev_potential(&bd, &sd);

      }
    }
  }
  return;    
}
  
}
