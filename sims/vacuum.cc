#include "vacuum.h"

namespace cosmo
{

void VacuumSim::init()
{
  t_init->start();

  // initialize base class
  simInit();

  cosmo_vacuum_db = input_db->getDatabase("VacuumSim");
    
  //  iodata->log("Running 'vacuum' type simulation.");

  pbox::pout<<"Running 'vacuum' type simulation."<<endl;
  t_init->stop();
}

/**
 * @brief      Set vacuum initial conditions
 *
 * @param[in]  map to BSSN fields
 * @param      initialized IOData
 */
void VacuumSim::setICs()
{
  iodata->log("Setting initial conditions (ICs).");

  gridding_algorithm->printClassData(tbox::plog);
  gridding_algorithm->makeCoarsestLevel();
  iodata->log("Finished setting ICs.");
}

void VacuumSim::initVacuumStep(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy))
{
  bssnSim->stepInit(hierarchy);
}

void VaccumSim::registerVariablesWithPlotter(
  appu::VisItDataWriter& visit_writer)
{
  std::vector<std::string> output_list =
    cosmo_vacuum_db->getStringVector("output_list");
  for (std::vector<std::string>::iterator it = output_list.begin() ;
       it != output_list.end(); ++it)
  {
    
    idx_t idx =
      variable_db->mapVariableAndContextToIndex(
        variable_db->getVariable(*it), variable_db->getContext("PREVIOUS"));
    
    if(idx < 0)
      TBOX_ERROR("Cannot find variable" <<*it<<" in output list!\n");
      
    visit_writer.registerPlotQuantity(
      *it,
      "SCALAR",
      idx,
      0,
      1.0,
      "CELL");
  }
}
  
void VacuumSim::initCoarsest(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  std::string ic_type = cosmo_vacuum_db->getString("ic_type");

  TBOX_ASSERT(ic_type);


  if(ic_type == "static_blackhole")
  {
    bssn_ic_static_blackhole(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
  }
  else
    TBOX_ERROR("Undefined IC type!\n");
  
  //already set _a = _p and _f = 0 

}

void VacuumSim::computeVectorWeights(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  TBOX_ASSERT(hierarchy);
  TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *hierarchy);

  int weight_id = weight_idx;
  int coarsest_ln = 0;
  int finest_ln = hierarchy->getFinestLevelNumber();


  int ln;
  for (ln = finest_ln; ln >= coarsest_ln; --ln) {

    /*
     * On every level, first assign cell volume to vector weight.
     */

    boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    for (hier::PatchLevel::iterator p(level->begin());
         p != level->end(); ++p) {
      const boost::shared_ptr<hier::Patch>& patch = *p;
      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      TBOX_ASSERT(patch_geometry);

      const double* dx = patch_geometry->getDx();
      double cell_vol = dx[0];
      if (d_dim > tbox::Dimension(1)) {
        cell_vol *= dx[1];
      }

      if (d_dim > tbox::Dimension(2)) {
        cell_vol *= dx[2];
      }

      boost::shared_ptr<pdat::CellData<double> > w(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_id)));
      TBOX_ASSERT(w);
      w->fillAll(cell_vol);
    }

    /*
     * On all but the finest level, assign 0 to vector
     * weight to cells covered by finer cells.
     */

    if (ln < finest_ln) {

      /*
       * First get the boxes that describe index space of the next finer
       * level and coarsen them to describe corresponding index space
       * at this level.
       */

      boost::shared_ptr<hier::PatchLevel> next_finer_level(
        hierarchy->getPatchLevel(ln + 1));
      hier::BoxContainer coarsened_boxes = next_finer_level->getBoxes();
      hier::IntVector coarsen_ratio(next_finer_level->getRatioToLevelZero());
      coarsen_ratio /= level->getRatioToLevelZero();
      coarsened_boxes.coarsen(coarsen_ratio);

      /*
       * Then set vector weight to 0 wherever there is
       * a nonempty intersection with the next finer level.
       * Note that all assignments are local.
       */

      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {

        const boost::shared_ptr<hier::Patch>& patch = *p;
        for (hier::BoxContainer::iterator i = coarsened_boxes.begin();
             i != coarsened_boxes.end(); ++i) {

          hier::Box intersection = *i * (patch->getBox());
          if (!intersection.empty()) {
            boost::shared_ptr<pdat::CellData<double> > w(
              BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(weight_id)));
            TBOX_ASSERT(w);
            w->fillAll(0.0, intersection);

          }  // assignment only in non-empty intersection
        }  // loop over coarsened boxes from finer level
      }  // loop over patches in level
    }  // all levels except finest
  }  // loop over levels
  
}

  
void VacuumSim::initializeLevelData(
   /*! Hierarchy to initialize */
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   /*! Level to initialize */
   const int ln,
   const double init_data_time,
   const bool can_be_refined,
   /*! Whether level is being introduced for the first time */
   const bool initial_time,
   /*! Level to copy data from */
   const boost::shared_ptr<hier::PatchLevel>& old_level,
   const bool allocate_data)
{
   NULL_USE(init_data_time);
   NULL_USE(can_be_refined);
   NULL_USE(initial_time);

   boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy(hierarchy);

   /*
    * Reference the level object with the given index from the hierarchy.
    */
   boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

   /*
    * If instructed, allocate all patch data on the level.
    * Allocate only persistent data.  Scratch data will
    * generally be allocated and deallocated as needed.
    */
   
   boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
       patch_hierarchy->getGridGeometry()));
   TBOX_ASSERT(grid_geometry_);
   geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

   if (allocate_data)
   {
     BSSN_APPLY_TO_FIELDS(RK4_ARRAY_ALLOC);
     level->allocatePatchData(weight_idx);
   }
   if(initial_time && ln == 0)
   {
     initCoarsest(patch_hierarchy);
   }
     
   /*
    * Refine solution data from coarser level and, if provided, old level.
    */
   
   xfer::RefineAlgorithm refiner;

   boost::shared_ptr<hier::RefineOperator> accurate_refine_op =
     bssnSim->space_refine_op;
     
   TBOX_ASSERT(accurate_refine_op);

   //registering refine variables
   BSSN_APPLY_TO_FIELDS(VAC_REGISTER_SPACE_REFINE_P);
      
   boost::shared_ptr<xfer::RefineSchedule> refine_schedule;

   level->getBoxLevel()->getMPI().Barrier();
   if (ln > 0) {
     /*
      * Include coarser levels in setting data
      */
     refine_schedule =
       refiner.createSchedule(
         level,
         old_level,
         ln - 1,
         hierarchy,
         NULL);
      
     else
     {
       /*
        * There is no coarser level, and source data comes only
        * from old_level, if any.
        */
       if (old_level)
       {
         refine_schedule =
           refiner.createSchedule(level,
                                  old_level,
                                  NULL);
       }
       else
         refine_schedule =
           refiner.createSchedule(level,
                                  level,
                                  NULL);
     }
     level->getBoxLevel()->getMPI().Barrier();
     

     if (refine_schedule)
     {
       refine_schedule->fillData(0.0);
         // It is null if this is the bottom level.
     }
     else
     {
       TBOX_ERROR(
         "Can not get refine schedule, check your code!\n");
     }

     math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);

     BSSN_APPLY_TO_FIELDS(BSSN_SWAP_PF);
     BSSN_APPLY_TO_FIELDS(BSSN_COPY_P_TO_A);
     BSSN_APPLY_TO_FIELDS(BSSN_SET_F_ZERO);

      
     if (0)
     {
       // begin debug code
       math::HierarchyCellDataOpsReal<double> hcellmath_debug(hierarchy);
       //hcellmath_debug.printData(d_phi_current, tbox::pout, false);
       // end debug code
     }
   }
   level->getBoxLevel()->getMPI().Barrier();
   /* Set vector weight. */
   computeVectorWeights(hierarchy);
}

void VacuumSim::applyGradientDetector(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy_,
   const int ln,
   const double error_data_time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation)
{
   NULL_USE(uses_richardson_extrapolation);
   NULL_USE(error_data_time);
   NULL_USE(initial_time);

   if (lstream) {
      *lstream
      << "VaccumSim("  << ")::applyGradientDetector"
      << std::endl;
   }
   hier::PatchHierarchy& hierarchy = *hierarchy_;
   boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
       hierarchy.getGridGeometry()));
   double max_der_norm = 0;
   hier::PatchLevel& level =
      (hier::PatchLevel &) * hierarchy.getPatchLevel(ln);
   size_t ntag = 0, ntotal = 0;
   double maxestimate = 0;
   for (hier::PatchLevel::iterator pi(level.begin());
        pi != level.end(); ++pi)
   {
      hier::Patch& patch = **pi;

      const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch.getPatchGeometry()));

      boost::shared_ptr<hier::PatchData> tag_data(
         patch.getPatchData(tag_index));
      ntotal += patch.getBox().numberCells().getProduct();
      if (!tag_data)
      {
         TBOX_ERROR(
            "Data index " << tag_index << " does not exist for patch.\n");
      }
      boost::shared_ptr<pdat::CellData<int> > tag_cell_data_(
         BOOST_CAST<pdat::CellData<int>, hier::PatchData>(tag_data));
      TBOX_ASSERT(tag_cell_data_);
      
      boost::shared_ptr<pdat::CellData<double>> K_data(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
          patch.getPatchData(bssnSim->DIFFK_p_idx)));

      boost::shared_ptr<pdat::CellData<double> > weight_(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
          patch.getPatchData(weight_idx)));

      arr_t K = pdat::ArrayDataAccess::access<DIM, real_t>(
        K_data->getArrayData())
      
      if (!K_data) {
         TBOX_ERROR("Data index " << bssnSim->K_p_idx
                                  << " does not exist for patch.\n");
      }
      pdat::CellData<idx_t>& tag_cell_data = *tag_cell_data_;
      pdat::CellData<real_t> & weight = *weight_;
          
      tag_cell_data.fill(0);
      
      hier::Box::iterator iend(patch.getBox().end());

      for (hier::Box::iterator i(patch.getBox().begin()); i != iend; ++i)
      {
         const pdat::CellIndex cell_index(*i);
         if(derivative_norm(
              cell_index(0),
              cell_index(1),
              cell_index(2),
              K) > d_adaption_threshold)
         {
           max_der_norm = derivative_norm(
              cell_index(0),
              cell_index(1),
              cell_index(2),
              K);
           tag_cell_data(cell_index) = 1;
           ++ntag;
         }
         
      }

   }

   tbox::plog << "Adaption threshold is " << d_adaption_threshold << "\n";
   tbox::plog << "Number of cells tagged on level " << ln << " is "
              << ntag << "/" << ntotal << "\n";
   tbox::plog << "Max norm is " << max_der_norm << "\n";
}
  
void VacuumSim::outputVacuumStep()
{
    // prepBSSNOutput();
    // io_bssn_fields_snapshot(iodata, step, bssnSim->fields);
    // io_bssn_fields_powerdump(iodata, step, bssnSim->fields, fourier);
    // io_bssn_dump_statistics(iodata, step, bssnSim->fields, bssnSim->frw);
    // io_bssn_constraint_violation(iodata, step, bssnSim);
}

void VacuumSim::runVacuumStep(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double from_t, double to_t))
{
  t_RK_steps->start();
    // Full RK step minus init()
  advanceLevel(hierarchy,
               0,
               from_t,
               to_t);

  t_RK_steps->stop();
}

double VaccumSim::getDt(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  return grid_geometry->getDx()[0] * dt_frac;
}

void VacuumSim::runStep(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  runCommonStepTasks(hierarcy);

  double dt = getDt(hierarchy);
  initVacuumStep(hierarchy);
  outputVacuumStep(hierarchy);
  runVacuumStep(hierarchy, cur_t, cur_t + dt);
  cur_t += dt;
}


void VaccumSim::RKEvolveLevel(
  const boost::shared_ptr<hier::PatchLevel> & level,
  double from_t,
  double to_t)
{
  xfer::RefineAlgorithm refiner;
  boost::shared_ptr<xfer::RefineSchedule> refine_schedule;

  bssnSim->registerSameLevelRefinerActive(refiner, space_refine_op);
  
  //BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REGISTER_SAME_LEVEL_REFINE_A,refiner,space_refine_op);

  refine_schedule = refiner.createSchedule(level, NULL);
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;

    //BSSN_APPLY_TO_FIELDS(BSSN_PDATA_ALL_INIT);
    //BSSN_APPLY_TO_FIELDS(BSSN_MDA_ACCESS_ALL_INIT);

    bssnSim->initPData(patch);
    bssnSim->initMDA(patch);
    
    bssnSim->RKEvolvePatch(patch);
    bssnSim->K1FinalizePatch(patch, to_t - from_t);  
  }

  level->getBoxLevel()->getMPI().Barrier();
  refine_schedule->fillData(to_t);

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    patch_geom = 
      BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch.getPatchGeometry());

    //BSSN_APPLY_TO_FIELDS(BSSN_PDATA_ALL_INIT);
    //BSSN_APPLY_TO_FIELDS(BSSN_MDA_ACCESS_ALL_INIT);

    bssnSim->initPData(patch);
    bssnSim->initMDA(patch);
    
    bssnSim->RKEvolvePatch(patch);
    bssnSim->K2FinalizePatch(patch, to_t - from_t);
  }

  
  level->getBoxLevel()->getMPI().Barrier();
  refine_schedule->fillData(to_t);

  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    patch_geom = 
      BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch.getPatchGeometry());

    //BSSN_APPLY_TO_FIELDS(BSSN_PDATA_ALL_INIT);
    //BSSN_APPLY_TO_FIELDS(BSSN_MDA_ACCESS_ALL_INIT);

    bssnSim->initPData(patch);
    bssnSim->initMDA(patch);
    
    bssnSim->RKEvolvePatch(patch);
    bssnSim->K3FinalizePatch(patch, to_t - from_t);
    
  }

  level->getBoxLevel()->getMPI().Barrier();
  refine_schedule->fillData(to_t);

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    patch_geom = 
      BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch.getPatchGeometry());

    //BSSN_APPLY_TO_FIELDS(BSSN_PDATA_ALL_INIT);
    //BSSN_APPLY_TO_FIELDS(BSSN_MDA_ACCESS_ALL_INIT);

    bssnSim->initPData(pathc);
    bssnSim->initMDA(patch);
    
    bssnSim->RKEvolvePatch(patch);
    bssnSim->K4FinalizePatch(patch, to_t - from_t);  
  }

  //BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REGISTER_SAME_LEVEL_REFINE_F,refiner,space_refine_op);

  bssnSim->registerSameLevelRefinerFinal(refiner,space_refine_op);
  
  refine_schedule = refiner.createSchedule(level, NULL);

  level->getBoxLevel()->getMPI().Barrier();
  refine_schedule->fillData(to_t);
}

void VacuumSim::advanceLevel(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  int ln,
  double from_t,
  double to_t)
{
  if( ln >= hierarchy->getNumberOfLevels())
    return;
  
  double dt = to_t - from_t;
  const boost::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  //RK advance interior(including innner ghost cells) of level
  RKEvolveLevel(level, from_t, to_t);

  setLevelTime(level, from_t, to_t);
  
  if(ln > 0)
  {
    boost::shared_ptr<xfer::PatchLevelBorderFillPattern> border_fill_pattern (
      new xfer::PatchLevelBorderFillPattern());

    xfer::RefineAlgorithm refiner;
    boost::shared_ptr<xfer::RefineSchedule> refine_schedule;

    bssnSim->registerInterLevelRefinerFinal(
      refiner,
      space_refine_op,
      time_refine_op);
    //BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REGISTER_INTER_LEVEL_REFINE_F,refiner,space_refine_op,time_refine_op);

    refine_schedule = refiner.createSchedule(
      border_fill_pattern,
      level,
      level,
      ln-1
      hierarchy,
      (cosmoPS->is_time_depent)?NULL:PS,
      TRUE,
      NULL);

    level->getBoxLevel()->getMPI().Barrier();
    refine_schedule->fillData(to_t);

  }
  
  advanceLevel(hierarchy, ln+1, from_t, from_t + (to_t - from_t)/2.0);

  advanceLevel(hierarchy, ln+1, from_t + (to_t - from_t)/2.0, to_t);


  if(ln < hierarchy->getNumberOfLevels() -1 )
  {
    xfer::CoarsenAlgorithm coarsener(dim);
  

    boost::shared_ptr<xfer::CoarsenSchedule> coarsen_schedule;

    //BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REGISTER_COARSEN_F, coarsener, coarsen_op);
    bssnSim->registerCoarsenFinal()
      
    coarsen_schedule = coarsener.createSchedule(level, hierarchy->getPatchLevel(ln+1));
    level->getBoxLevel()->getMPI().Barrier();
    coarsen_schedule->coarsenData();

    xfer::RefineAlgorithm post_refiner;

    boost::shared_ptr<xfer::RefineSchedule> refine_schedule;

    
    //BSSN_APPLY_TO_FIELDS_ARGS(BSSN_REGISTER_SAME_LEVEL_REFINE_F,post_refiner,space_refine_op);

    bssnSim->registerSameLevelRefiner(post_refiner, space_refine_op);
    
    refine_schedule = post_refiner.createSchedule(level, NULL);

    level->getBoxLevel()->getMPI().Barrier();
    refine_schedule->fillData(to_t);
  }

  // swap _p and _f patches, so all recent data is stored
  // in _p arrays rather than _f arrays 

  math::HierarchyCellDataOpsReal<real_t> hcellmath(hierarchy,ln,ln);

  bssnSim->swapPF(hcellmath);
  bssnSim->copyPToA(hcellmath);
  bssnSim->setFZero(hcellmath);
  
}
  


} /* namespace cosmo */
