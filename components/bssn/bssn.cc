#include "bssn.h"
#include "../../cosmo_includes.h"
#include "../../utils/math.h"


using namespace SAMRAI;

namespace cosmo
{

/**
 * @brief Constructor for BSSN class
 * @details Allocate memory for fields, add fields to map,
 * create reference FRW integrator, and call BSSN::init.
 */
BSSN::BSSN(
  const tbox::Dimension& dim_in,
  boost::shared_ptr<tbox::Database> database_in,
  std::ostream* l_stream_in,
  real_t KO_damping_coefficient_in):
  lstream(l_stream_in),
  cosmo_bssn_db(database_in),
  dim(dim_in),
  KO_damping_coefficient(KO_damping_coefficient_in),
  gaugeHandler(new BSSNGaugeHandler(cosmo_bssn_db)),
  g_eta(cosmo_bssn_db->getDoubleWithDefault("g_eta",1.0))
{
  
  BSSN_APPLY_TO_FIELDS(VAR_INIT);
  BSSN_APPLY_TO_SOURCES(VAR_INIT);
  BSSN_APPLY_TO_GEN1_EXTRAS(VAR_INIT);

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


  
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REG_TO_CONTEXT, context_scratch, s, STENCIL_ORDER);
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REG_TO_CONTEXT, context_previous, p, STENCIL_ORDER);
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REG_TO_CONTEXT, context_active, a, STENCIL_ORDER);
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REG_TO_CONTEXT, context_k1, k1, STENCIL_ORDER);
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REG_TO_CONTEXT, context_k2, k2, STENCIL_ORDER);
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REG_TO_CONTEXT, context_k3, k3, STENCIL_ORDER);
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REG_TO_CONTEXT, context_k4, k4, STENCIL_ORDER);


  BSSN_APPLY_TO_SOURCES_ARGS(BSSN_REG_TO_CONTEXT, context_active, a, STENCIL_ORDER);

  BSSN_APPLY_TO_GEN1_EXTRAS_ARGS(BSSN_REG_TO_CONTEXT, context_active, a, STENCIL_ORDER);

  init();
  
}

BSSN::~BSSN()
{
  // BSSN_APPLY_TO_FIELDS(RK4_ARRAY_DELETE)
  // BSSN_APPLY_TO_SOURCES(GEN1_ARRAY_DELETE)
  // BSSN_APPLY_TO_GEN1_EXTRAS(GEN1_ARRAY_DELETE)
}

/**
 * @brief Initialize fields in BSSN class to defaults
 * @details BSSN fields initialized to a flat (difference) metric with zero source;
 * thus all fields in all registers are zeroed. Reference integrator unaffected.
 */
void BSSN::init()
{
  
}


/**
 * @brief Call RK4Register class step initialization; normalize Aij and DIFFgammaIJ fields
 * @details See RK4Register::stepInit() method.
 */
void BSSN::stepInit(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  
  //BSSN_RK_INITIALIZE; // macro calls stepInit for all fields

    boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
       hierarchy->getGridGeometry()));
   TBOX_ASSERT(grid_geometry_);
   geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  const double * domain_lower = &grid_geometry.getXLower()[0];
  const double * domain_upper = &grid_geometry.getXUpper()[0];

  for(int i = 0 ; i < DIM; i++)
    L[i] = domain_upper[i] - domain_lower[i];

  
// # if NORMALIZE_GAMMAIJ_AIJ
//     set_DIFFgamma_Aij_norm(); // norms _a register
// # endif
}

void BSSN::RKEvolvePatchBD(
  const boost::shared_ptr<hier::Patch> & patch,
  real_t dt)
{
  
  boost::shared_ptr<hier::PatchGeometry> geom (patch->getPatchGeometry());

  idx_t codim = 1;
  const std::vector<hier::BoundaryBox> & codim1_boxes =
    geom->getCodimensionBoundaries(codim);

  const idx_t n_codim1_boxes = static_cast<idx_t>(codim1_boxes.size());

  if(n_codim1_boxes == 0) return;

  initPData(patch);
  initMDA(patch);

  const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom( 
    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

    
  //initialize dx for each patch
  //const real_t * dx = &(patch_geom->getDx())[0];

  const hier::Box& patch_box = patch->getBox();

  //hier::Box & boundary_fill_box;

  for(int i = 0 ; i < n_codim1_boxes; i++)
  {
    hier::Box boundary_fill_box =
      geom->getBoundaryFillBox(
        codim1_boxes[i], patch_box, DIFFchi_a_pdata->getGhostCellWidth());

    if(boundary_fill_box.empty()) continue;


    idx_t l_idx = codim1_boxes[i].getLocationIndex();
    
    BSSNData bd = {0};

    //tbox::pout<<boundary_fill_box<<"%%%%%%\n";
    boundary_fill_box.shift(
      (hier::Box::dir_t)l_idx/2,
      (l_idx%2)?(-STENCIL_ORDER_WIDTH):STENCIL_ORDER_WIDTH);


    boundary_fill_box *= patch_box;
    //tbox::pout<<boundary_fill_box<<"*******\n";
    //    throw(-1);
    //initialize dx for each patch
    const real_t * dx = &(patch_geom->getDx())[0];

    const idx_t * lower = &boundary_fill_box.lower()[0];
    const idx_t * upper = &boundary_fill_box.upper()[0];

    
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          set_bd_values_bd(i, j, k, &bd, dx);
          BSSN_RK_EVOLVE_BD;
        }
      }
    }
    //    tbox::pout<<boundary_fill_box<<"\n";
  }
  /************************updating codim = 2 boundaries****************/
  codim = 2;

  const std::vector<hier::BoundaryBox> & codim2_boxes =
    geom->getCodimensionBoundaries(codim);

  const idx_t n_codim2_boxes = static_cast<idx_t>(codim2_boxes.size());




  for(int i = 0 ; i < n_codim2_boxes; i++)
  {
    hier::Box  boundary_fill_box =
      geom->getBoundaryFillBox(
        codim2_boxes[i], patch_box, DIFFchi_a_pdata->getGhostCellWidth());
    
    if(boundary_fill_box.empty()) continue;  


    idx_t l_idx = codim2_boxes[i].getLocationIndex();
    
    BSSNData bd = {0};

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

    //tbox::pout<<boundary_fill_box<<"\n";
    boundary_fill_box.shift(hier::IntVector(shift_vec));

    boundary_fill_box *= patch_box;
    //tbox::pout<<boundary_fill_box<<"\n";

    //initialize dx for each patch
    const real_t * dx = &(patch_geom->getDx())[0];

    const idx_t * lower = &boundary_fill_box.lower()[0];
    const idx_t * upper = &boundary_fill_box.upper()[0];

    
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
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




  for(int i = 0 ; i < n_codim3_boxes; i++)
  {
    hier::Box boundary_fill_box =
      geom->getBoundaryFillBox(
        codim3_boxes[i], patch_box, DIFFchi_a_pdata->getGhostCellWidth());

    if(boundary_fill_box.empty()) continue;

    idx_t l_idx = codim3_boxes[i].getLocationIndex();

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

    
    BSSNData bd = {0};


    //initialize dx for each patch
    const real_t * dx = &(patch_geom->getDx())[0];
  
    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          set_bd_values_bd(i, j, k, &bd, dx);
          BSSN_RK_EVOLVE_BD;
        }
      }
    }
  }
  
}

#if USE_CCZ4
void BSSN::initZ(
  const boost::shared_ptr<hier::PatchLevel> & level)
{
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    const hier::Box& box = DIFFchi_a_pdata->getBox();
    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];
    BSSNData bd = {0};

    initPData(patch);
    initMDA(patch);
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
          set_bd_values_for_extra_fields(i,j,k,&bd,dx);
          Z1_a(i,j,k) =
            (bd.gamma11*(bd.Gamma1 - bd.Gammad1)
             + bd.gamma12*(bd.Gamma2 - bd.Gammad2)
             + bd.gamma13*(bd.Gamma3 - bd.Gammad3))/2.0;

          Z2_a(i,j,k) =
            (bd.gamma21*(bd.Gamma1 - bd.Gammad1)
             + bd.gamma22*(bd.Gamma2 - bd.Gammad2)
             + bd.gamma23*(bd.Gamma3 - bd.Gammad3))/2.0;

          Z3_a(i,j,k) =
            (bd.gamma31*(bd.Gamma1 - bd.Gammad1)
             + bd.gamma32*(bd.Gamma2 - bd.Gammad2)
             + bd.gamma33*(bd.Gamma3 - bd.Gammad3))/2.0;
        }
      }
    }
    
  }

}
void BSSN::initZ(
  const boost::shared_ptr<hier::Patch> & patch)
{
  
  const hier::Box& box = DIFFchi_a_pdata->getGhostBox();
  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];
  BSSNData bd = {0};

  initPData(patch);
  initMDA(patch);
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
        set_bd_values_for_extra_fields(i,j,k,&bd,dx);
        Z1_a(i,j,k) =
          (bd.gamma11*(bd.Gamma1 - bd.Gammad1)
           + bd.gamma12*(bd.Gamma2 - bd.Gammad2)
           + bd.gamma13*(bd.Gamma3 - bd.Gammad3))/2.0;

        Z2_a(i,j,k) =
          (bd.gamma21*(bd.Gamma1 - bd.Gammad1)
           + bd.gamma22*(bd.Gamma2 - bd.Gammad2)
           + bd.gamma23*(bd.Gamma3 - bd.Gammad3))/2.0;

        Z3_a(i,j,k) =
          (bd.gamma31*(bd.Gamma1 - bd.Gammad1)
           + bd.gamma32*(bd.Gamma2 - bd.Gammad2)
           + bd.gamma33*(bd.Gamma3 - bd.Gammad3))/2.0;
      }
    }
  }
  

}

#endif
/**
 * @brief Call BSSN::RKEvolvePt for all points
 */
void BSSN::RKEvolvePatch(
  const boost::shared_ptr<hier::Patch> & patch, real_t dt)
{
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
        set_bd_values(i, j, k, &bd, dx);
        BSSN_RK_EVOLVE_PT;
      }
    }
  }
}


void BSSN::prepairForK1(
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
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(DIFFchi_a_idx)->getTime()) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(BSSN_INIT_R_K1);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(DIFFchi_p_idx)->getTime())
       - (patch->getPatchData(DIFFchi_a_idx)->getTime()
          - patch->getPatchData(DIFFchi_p_idx)->getTime())) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(BSSN_INIT_L_K1);
          }
        }
      }
    }
    else
      TBOX_ERROR("Current level locates neigher L or R branch of its father, check your code!");
    
  }
}

void BSSN::prepairForK2(
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
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(DIFFchi_a_idx)->getTime()) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(BSSN_INIT_R_K2);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(DIFFchi_p_idx)->getTime())
       - (patch->getPatchData(DIFFchi_a_idx)->getTime()
          - patch->getPatchData(DIFFchi_p_idx)->getTime())) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(BSSN_INIT_L_K2);
          }
        }
      }
    }
  }
}

void BSSN::prepairForK3(
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

    BSSNData bd = {0};

    // right branch of the tree
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(DIFFchi_a_idx)->getTime()) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(BSSN_INIT_R_K3);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(DIFFchi_p_idx)->getTime())
       - (patch->getPatchData(DIFFchi_a_idx)->getTime()
          - patch->getPatchData(DIFFchi_p_idx)->getTime())) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(BSSN_INIT_L_K3);
          }
        }
      }
    }
  }
}

//final step
void BSSN::prepairForK4(
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

    
    
    BSSNData bd = {0};


    // right branch of the tree
    if(tbox::MathUtilities<real_t>::Abs(to_t - patch->getPatchData(DIFFchi_a_idx)->getTime()) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(BSSN_INIT_R_K4);
          }
        }
      }
    }

    // left branch of the tree
    else if(tbox::MathUtilities<real_t>::Abs(2.0*(to_t - patch->getPatchData(DIFFchi_p_idx)->getTime())
       - (patch->getPatchData(DIFFchi_a_idx)->getTime()
          - patch->getPatchData(DIFFchi_p_idx)->getTime())) < EPS)
    {
      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            BSSN_APPLY_TO_FIELDS(BSSN_INIT_L_K4);
          }
        }
      }
    }
  }
}  
  
void BSSN::registerRKRefiner(
  xfer::RefineAlgorithm& refiner,
  boost::shared_ptr<hier::RefineOperator> &space_refine_op)
{
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REGISTER_RK_REFINER, refiner, space_refine_op);
}

void BSSN::registerSameLevelRefinerActive(
  xfer::RefineAlgorithm& refiner,
  boost::shared_ptr<hier::RefineOperator> &space_refine_op)
{
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REGISTER_SAME_LEVEL_REFINE_A,refiner,space_refine_op);
}


// void BSSN::registerInterLevelRefinerActive(
//   xfer::RefineAlgorithm& refiner,
//   boost::shared_ptr<hier::RefineOperator> &space_refine_op
//   boost::shared_ptr<hier::TimeInterpolateOperator> & time_refine_op )
// {
//   BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REGISTER_INTER_LEVEL_REFINE_A,refiner,space_refine_op,time_refine_op);
// }

// void BSSN::registerInterLevelRefinerFinal(
//   xfer::RefineAlgorithm& refiner,
//   boost::shared_ptr<hier::RefineOperator>& space_refine_op)
//   boost::shared_ptr<hier::TimeInterpolateOperator> & time_refine_op )
// {
//   BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REGISTER_INTER_LEVEL_REFINE_F,refiner,space_refine_op,time_refine_op);
// }


void BSSN::registerCoarsenActive(
  xfer::CoarsenAlgorithm& coarsener,
  boost::shared_ptr<hier::CoarsenOperator>& coarsen_op)
{
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REGISTER_COARSEN_A, coarsener, coarsen_op);
}






void BSSN::copyAToP(
  math::HierarchyCellDataOpsReal<real_t> & hcellmath)
{
  BSSN_APPLY_TO_FIELDS(BSSN_COPY_A_TO_P);
}




void BSSN::initPData(
  const boost::shared_ptr<hier::Patch> & patch)
{
  BSSN_APPLY_TO_FIELDS(BSSN_PDATA_ALL_INIT);
  BSSN_APPLY_TO_SOURCES_ARGS(BSSN_PDATA_INIT, a);
  BSSN_APPLY_TO_GEN1_EXTRAS_ARGS(BSSN_PDATA_INIT, a);
}

void BSSN::initMDA(
  const boost::shared_ptr<hier::Patch> & patch)
{
  BSSN_APPLY_TO_FIELDS(BSSN_MDA_ACCESS_ALL_INIT);

  BSSN_APPLY_TO_SOURCES_ARGS(BSSN_MDA_ACCESS_INIT, a);
  
  BSSN_APPLY_TO_GEN1_EXTRAS_ARGS(BSSN_MDA_ACCESS_INIT, a);
}


    
void BSSN::setLevelTime(
  const boost::shared_ptr<hier::PatchLevel> & level,
  double from_t, double to_t)
{
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    BSSN_APPLY_TO_FIELDS_ARGS(BSSN_SET_PATCH_TIME, from_t, to_t);
  }
}
  
void BSSN::K1FinalizePatch(
  const boost::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);
  
  const hier::Box& box = DIFFchi_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        BSSN_FINALIZE_K(1);
        // if(i == 60 && j == 56 && k == 63 && patch->getPatchLevelNumber() == 2)
        //   tbox::pout<<DIFFchi_p(i,j,k)<<"##"<<patch->getGlobalId()<<"\n";
      }
    }
  }
}

void BSSN::K2FinalizePatch(
  const boost::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  const hier::Box& box = DIFFchi_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
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
  const boost::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  const hier::Box& box = DIFFchi_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
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


void BSSN::K4FinalizePatch(
  const boost::shared_ptr<hier::Patch> & patch)
{
  initPData(patch);
  initMDA(patch);

  //Only evolv inner grids for K4
  const hier::Box& box = DIFFchi_a_pdata->getGhostBox();

  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
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
}

void BSSN::set_bd_values_for_extra_fields(
  idx_t i, idx_t j, idx_t k, BSSNData *bd, const real_t dx[])
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
  
  // draw data from cache
  //set_local_vals(bd);
  BSSN_APPLY_TO_FIELDS(RK4_SET_LOCAL_VALUES);
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
  
  calculate_dgamma(bd, dx);
  calculate_conformal_christoffels(bd, dx);
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
# if USE_CCZ4
  calculate_dtheta(bd,dx);
# endif
# if USE_BSSN_SHIFT
  calculate_dbeta(bd,dx);
# endif

  #if USE_EXPANSION
  calculate_dexpN(bd,dx);
  #endif
  // Christoffels depend on metric & derivs.
  calculate_conformal_christoffels(bd,dx);
  // DDw depend on christoffels, metric, and derivs
  calculateDDphi(bd,dx);
  calculateDDalphaTF(bd,dx);
  // Ricci depends on DDphi
  calculateRicciTF(bd,dx);

  calculateDZ(bd,dx);
  // Hamiltonian constraint
  //bd->H = hamiltonianConstraintCalc(bd,dx);
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

  bd->d1phi = - 0.5 * bd->d1chi / bd->chi;
  bd->d2phi = - 0.5 * bd->d2chi / bd->chi;
  bd->d3phi = - 0.5 * bd->d3chi / bd->chi;

  bd->d1d1phi = 0.5*(bd->d1chi*bd->d1chi/pw2(bd->chi)-bd->d1d1chi/bd->chi); 
  bd->d1d2phi = 0.5*(bd->d1chi*bd->d2chi/pw2(bd->chi)-bd->d1d2chi/bd->chi);
  bd->d1d3phi = 0.5*(bd->d1chi*bd->d3chi/pw2(bd->chi)-bd->d1d3chi/bd->chi);
  bd->d2d2phi = 0.5*(bd->d2chi*bd->d2chi/pw2(bd->chi)-bd->d2d2chi/bd->chi);
  bd->d2d3phi = 0.5*(bd->d2chi*bd->d3chi/pw2(bd->chi)-bd->d2d3chi/bd->chi);
  bd->d3d3phi = 0.5*(bd->d3chi*bd->d3chi/pw2(bd->chi)-bd->d3d3chi/bd->chi); 
  
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

#if USE_CCZ4
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

#if USE_EXPANSION
void BSSN::calculate_dexpN(BSSNData *bd, const real_t dx[])
{
  bd->d1expN = upwind_derivative(bd->i, bd->j, bd->k, 1, expN_a, dx, bd->beta1);
  bd->d2expN = upwind_derivative(bd->i, bd->j, bd->k, 2, expN_a, dx, bd->beta2);
  bd->d3expN = upwind_derivative(bd->i, bd->j, bd->k, 3, expN_a, dx, bd->beta3);
}
#endif


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

void BSSN::calculateDDphi(BSSNData *bd, const real_t dx[])
{
  // double covariant derivatives, using unitary metric
  bd->D1D1phi = bd->d1d1phi - (bd->G111*bd->d1phi + bd->G211*bd->d2phi + bd->G311*bd->d3phi);
  bd->D2D2phi = bd->d2d2phi - (bd->G122*bd->d1phi + bd->G222*bd->d2phi + bd->G322*bd->d3phi);
  bd->D3D3phi = bd->d3d3phi - (bd->G133*bd->d1phi + bd->G233*bd->d2phi + bd->G333*bd->d3phi);

  bd->D1D2phi = bd->d1d2phi - (bd->G112*bd->d1phi + bd->G212*bd->d2phi + bd->G312*bd->d3phi);
  bd->D1D3phi = bd->d1d3phi - (bd->G113*bd->d1phi + bd->G213*bd->d2phi + bd->G313*bd->d3phi);
  bd->D2D3phi = bd->d2d3phi - (bd->G123*bd->d1phi + bd->G223*bd->d2phi + bd->G323*bd->d3phi);  
}

void BSSN::calculateDZ(BSSNData *bd, const real_t dx[])
{
  //
  // bd->D1Z1 = derivative(bd->i, bd->j, bd->k, 1, Z1_a, dx)
  //   - (bd->G111 * bd->Z1 + bd->G211 * bd->Z2 + bd->G311 * bd->Z3);
  // bd->D2Z2 = derivative(bd->i, bd->j, bd->k, 2, Z2_a, dx)
  //   - (bd->G122 * bd->Z1 + bd->G222 * bd->Z2 + bd->G322 * bd->Z3);
  // bd->D3Z3 = derivative(bd->i, bd->j, bd->k, 3, Z3_a, dx)
  //   - (bd->G133 * bd->Z1 + bd->G233 * bd->Z2 + bd->G333 * bd->Z3);
  
  // bd->D1Z2 = derivative(bd->i, bd->j, bd->k, 1, Z2_a, dx)
  //   - (bd->G112 * bd->Z1 + bd->G212 * bd->Z2 + bd->G312 * bd->Z3);
  // bd->D1Z3 = derivative(bd->i, bd->j, bd->k, 1, Z3_a, dx)
  //   - (bd->G113 * bd->Z1 + bd->G213 * bd->Z2 + bd->G313 * bd->Z3);
  // bd->D2Z3 = derivative(bd->i, bd->j, bd->k, 2, Z3_a, dx)
  //   - (bd->G123 * bd->Z1 + bd->G223 * bd->Z2 + bd->G323 * bd->Z3);

  // bd->D2Z1 = derivative(bd->i, bd->j, bd->k, 2, Z1_a, dx)
  //   - (bd->G121 * bd->Z1 + bd->G221 * bd->Z2 + bd->G321 * bd->Z3);
  // bd->D3Z1 = derivative(bd->i, bd->j, bd->k, 3, Z1_a, dx)
  //   - (bd->G131 * bd->Z1 + bd->G231 * bd->Z2 + bd->G331 * bd->Z3);
  // bd->D3Z2 = derivative(bd->i, bd->j, bd->k, 3, Z2_a, dx)
  //   - (bd->G132 * bd->Z1 + bd->G232 * bd->Z2 + bd->G332 * bd->Z3);

  bd->D1Z1 = bd->D1Z2 =bd->D1Z3=bd->D2Z1=bd->D2Z2=bd->D2Z3=bd->D3Z1=bd->D3Z2=bd->D3Z3=0;
  bd->DZTR = bd->gammai11 * bd->D1Z1 + bd->gammai22 * bd->D2Z2 + bd->gamma33 * bd->D3Z3
    + (bd->gammai12 * bd->D1Z2 + bd->gammai13 * bd->D1Z3 + bd->gamma23 * bd->D2Z3)
    + (bd->gammai21 * bd->D2Z1 + bd->gammai31 * bd->D3Z1 + bd->gamma32 * bd->D3Z2);

  bd->D1Z1TF = bd->D1Z1 - (1.0/3.0) * bd->gamma11 * bd->DZTR;
  bd->D2Z2TF = bd->D2Z2 - (1.0/3.0) * bd->gamma22 * bd->DZTR;
  bd->D3Z3TF = bd->D3Z3 - (1.0/3.0) * bd->gamma33 * bd->DZTR;

  bd->D1Z2TF = bd->D1Z2 - (1.0/3.0) * bd->gamma12 * bd->DZTR;
  bd->D1Z3TF = bd->D1Z3 - (1.0/3.0) * bd->gamma13 * bd->DZTR;
  bd->D2Z3TF = bd->D2Z3 - (1.0/3.0) * bd->gamma23 * bd->DZTR;

  bd->D2Z1TF = bd->D2Z1 - (1.0/3.0) * bd->gamma21 * bd->DZTR;
  bd->D3Z1TF = bd->D3Z1 - (1.0/3.0) * bd->gamma31 * bd->DZTR;
  bd->D3Z2TF = bd->D3Z2 - (1.0/3.0) * bd->gamma32 * bd->DZTR;

  bd->DZTR *= pw2(bd->chi);
}

void BSSN::calculateDDalphaTF(BSSNData *bd, const real_t dx[])
{
  // double covariant derivatives - use non-unitary metric - extra pieces that depend on phi!
  // the gammaIldlphi are needed for the BSSN_CALCULATE_DIDJALPHA macro
  real_t gammai1ldlphi = bd->gammai11*bd->d1phi + bd->gammai12*bd->d2phi + bd->gammai13*bd->d3phi;
  real_t gammai2ldlphi = bd->gammai21*bd->d1phi + bd->gammai22*bd->d2phi + bd->gammai23*bd->d3phi;
  real_t gammai3ldlphi = bd->gammai31*bd->d1phi + bd->gammai32*bd->d2phi + bd->gammai33*bd->d3phi;
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
    bd->gammai11*(bd->D1D1phi + 2.0*bd->d1phi*bd->d1phi)
    + bd->gammai22*(bd->D2D2phi + 2.0*bd->d2phi*bd->d2phi)
    + bd->gammai33*(bd->D3D3phi + 2.0*bd->d3phi*bd->d3phi)
    + 2.0*(
      bd->gammai12*(bd->D1D2phi + 2.0*bd->d1phi*bd->d2phi)
      + bd->gammai13*(bd->D1D3phi + 2.0*bd->d1phi*bd->d3phi)
      + bd->gammai23*(bd->D2D3phi + 2.0*bd->d2phi*bd->d3phi)
    )
  );

  bd->ricci11 = bd->Uricci11 - 2.0*( bd->D1D1phi - 2.0*bd->d1phi*bd->d1phi + bd->gamma11*(expression) );
  bd->ricci12 = bd->Uricci12 - 2.0*( bd->D1D2phi - 2.0*bd->d1phi*bd->d2phi + bd->gamma12*(expression) );
  bd->ricci13 = bd->Uricci13 - 2.0*( bd->D1D3phi - 2.0*bd->d1phi*bd->d3phi + bd->gamma13*(expression) );
  bd->ricci22 = bd->Uricci22 - 2.0*( bd->D2D2phi - 2.0*bd->d2phi*bd->d2phi + bd->gamma22*(expression) );
  bd->ricci23 = bd->Uricci23 - 2.0*( bd->D2D3phi - 2.0*bd->d2phi*bd->d3phi + bd->gamma23*(expression) );
  bd->ricci33 = bd->Uricci33 - 2.0*( bd->D3D3phi - 2.0*bd->d3phi*bd->d3phi + bd->gamma33*(expression) );
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
          1.0/3.0*(bd->K)*(bd->K)
#       endif

#       if !(EXCLUDE_SECOND_ORDER_SMALL)
          + bd->AijAij
#       endif
          + 2.0 * bd->DZTR
          - 2.0 * bd->theta * bd->K
    )
    + 4.0*PI*bd->alpha*(bd->r + bd->S)
    + upwind_derivative(bd->i, bd->j, bd->k, 1, DIFFK_a, dx, bd->beta1)
    + upwind_derivative(bd->i, bd->j, bd->k, 2, DIFFK_a, dx, bd->beta2)
    + upwind_derivative(bd->i, bd->j, bd->k, 3, DIFFK_a, dx, bd->beta3)
    - 3.0*bd->alpha*Z4c_K1_DAMPING_AMPLITUDE*(1.0 + Z4c_K2_DAMPING_AMPLITUDE)*bd->theta
    - KO_dissipation_Q(bd->i, bd->j, bd->k, DIFFK_a, dx, KO_damping_coefficient)
  );
}

//flag
real_t BSSN::ev_DIFFchi(BSSNData *bd, const real_t dx[])
{
  return (
    1.0/3.0*(
      bd->alpha* bd->chi * bd->K
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
  #if USE_CCZ4
  return (
    0.5*bd->alpha*(
      bd->ricci + 2.0/3.0*pw2( bd->K) - 2.0 * bd->theta * bd->K + 2.0 * bd->DZTR
      - bd->AijAij - 16.0*PI* bd->r)
    #if USE_BSSN_SHIFT
    + upwind_derivative(bd->i, bd->j, bd->k, 1, theta_a, dx, bd->beta1)
    + upwind_derivative(bd->i, bd->j, bd->k, 2, theta_a, dx, bd->beta2)
    + upwind_derivative(bd->i, bd->j, bd->k, 3, theta_a, dx, bd->beta3)
    #endif
    - bd->Z1 * bd->d1a - bd->Z2 * bd->d2a - bd->Z3 * bd->d3a
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
    -bd->alpha * bd->K/3.0;
}
#endif


#if USE_GAMMA_DRIVER
real_t BSSN::ev_auxB1(BSSNData *bd, const real_t dx[])
{
  return 0.75*ev_Gamma1(bd,dx)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 1, auxB1_a, dx, bd->beta1)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 2, auxB1_a, dx, bd->beta2)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 3, auxB1_a, dx, bd->beta3)
    - g_eta * bd->auxB1;
}

real_t BSSN::ev_auxB2(BSSNData *bd, const real_t dx[])
{
  return 0.75*ev_Gamma2(bd,dx)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 1, auxB2_a, dx, bd->beta1)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 2, auxB2_a, dx, bd->beta2)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 3, auxB2_a, dx, bd->beta3)
    - g_eta * bd->auxB2;
}

real_t BSSN::ev_auxB3(BSSNData *bd, const real_t dx[])
{
  return 0.75*ev_Gamma3(bd,dx)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 1, auxB3_a, dx, bd->beta1)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 2, auxB3_a, dx, bd->beta2)
    //+ upwind_derivative(bd->i, bd->j, bd->k, 3, auxB3_a, dx, bd->beta3)
    - g_eta * bd->auxB3;
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
            + bd->DIFFK           );
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
#if USE_CCZ4  
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


void BSSN::output_max_H_constaint(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t weight_idx)
{
  double max_H=0, max_H_scaled = 0;
  double max_M=0, max_M_scaled = 0;
  idx_t mp[3] = {0}, hp[3] = {0};  
  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln ++)
  {
    boost::shared_ptr <hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

    boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
       hierarchy->getGridGeometry()));
   TBOX_ASSERT(grid_geometry_);
   geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

    
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
    {
      const boost::shared_ptr<hier::Patch> & patch = *pit;

      const hier::Box& box = patch->getBox();

      const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      initPData(patch);

      initMDA(patch);


      boost::shared_ptr<pdat::CellData<double> > weight(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_idx)));
      

      arr_t weight_array =
        pdat::ArrayDataAccess::access<DIM, double>(
          weight->getArrayData());


      
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];
      
      const double *dx = &patch_geom->getDx()[0];

      
      BSSNData bd = {0};


      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            set_bd_values(i,j,k,&bd,dx);
            if(weight_array(i,j,k) > 0)
            {
              if(tbox::MathUtilities<double>::Abs(
                  hamiltonianConstraintCalc(&bd, dx)) > max_H)
              {
                hp[0] = i, hp[1] = j, hp[2] = k;
              }
              max_H = tbox::MathUtilities<double>::Max(
                max_H, tbox::MathUtilities<double>::Abs(
                  hamiltonianConstraintCalc(&bd, dx)));
               max_H_scaled = tbox::MathUtilities<double>::Max(
                max_H_scaled, tbox::MathUtilities<double>::Abs(
                  hamiltonianConstraintCalc(&bd, dx)/
                  hamiltonianConstraintScale(&bd,dx)));

               if(tbox::MathUtilities<double>::Abs(
                  momentumConstraintCalc(&bd, dx)) > max_M)
               {
                 mp[0]=i, mp[1] =j , mp[2] = k;
               }

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
            <<max_H<<"/"<<max_H_scaled<<"\n at position "<<hp[0]<<" "<<hp[1]<<" "<<hp[2]<<"\n";
  tbox::pout<<"Max Momentum constraint is "
            <<max_M<<"/"<<max_M_scaled<<"\n at position "<<mp[0]<<" "<<mp[1]<<" "<<mp[2]<<"\n";
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
      bd->ricci + 2.0/3.0*pw2(bd->K ) - bd->AijAij - 16.0*PI*bd->r
    );
}

real_t BSSN::hamiltonianConstraintScale(BSSNData *bd, const real_t dx[])
{

  // sqrt sum of sq. of terms for appx. mag / scale
  return 
    pow(bd->chi, -2.5)/8.0*sqrt( pw2(bd->ricci) + pw2(bd->AijAij)
          + pw2(2.0/3.0*pw2(bd->K ))
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
