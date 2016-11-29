#include "vacuum.h"

namespace cosmo
{

void VacuumSim::init()
{
  t_init->start();

  // initialize base class
  simInit();

  
    
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

  // meh
  // TODO: Set vacuum ICs (eg, AwA test)

  iodata->log("Finished setting ICs.");
}

void VacuumSim::initVacuumStep()
{
  bssnSim->stepInit();
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
      d_phi_previous,
      0,
      1.0,
      "CELL");

  }
}
  
void VacuumSim::initCoarsest(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  boost::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(0));

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
    
    boost::shared_ptr<pdat::CellData<double> > pi_previous(
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
        patch->getPatchData(d_pi_previous)));
    boost::shared_ptr<pdat::CellData<double> > phi_previous(
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
        patch->getPatchData(d_phi_previous)));


    MDA_Access<double, 3, MDA_OrderColMajor<3> > pi =
      pdat::ArrayDataAccess::access<3, double>(
        pi_previous->getArrayData());
    MDA_Access<double, 3, MDA_OrderColMajor<3> > phi =
      pdat::ArrayDataAccess::access<3, double>(
        phi_previous->getArrayData());

    const double * domain_lower = &grid_geometry.getXLower()[0];
    const double * domain_upper = &grid_geometry.getXUpper()[0];

    const double *dx = &grid_geometry.getDx()[0];

      
    double L[3];

    for(int i = 0 ; i < 3; i++)
      L[i] = domain_upper[i] - domain_lower[i];
    
    d_omega = sqrt( PW2(2.0 * PI * nx / L[0])
                    + PW2(2.0 * PI * ny / L[1])
                    + PW2(2.0 * PI * nz / L[2]));

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          phi(i,j,k) = sin(2.0 * PI * (double) i * dx[0] * nx / L[0] +
                           2.0 * PI * (double) j * dx[1] * ny / L[1] +
                           2.0 * PI * (double) k * dx[2] * nz / L[2]);
          
          pi(i,j,k) = d_omega *
            cos(2.0 * PI * (double) i * dx[0] * nx / L[0] +
                           2.0 * PI * (double) j * dx[1] * ny / L[1] +
                           2.0 * PI * (double) k * dx[2] * nz / L[2]);
        }
      }
    }

        
  }
  

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
   {
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
      }
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
     math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);
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
      
      if (0) {
         // begin debug code
         math::HierarchyCellDataOpsReal<double> hcellmath_debug(hierarchy);
         hcellmath_debug.printData(d_phi_current, tbox::pout, false);
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
          patch.getPatchData(bssnSim->K_p_idx)));

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

void VacuumSim::runVacuumStep()
{
  t_RK_steps->start();
    // Full RK step minus init()
  bssnSim->step();

  t_RK_steps->stop();
}

void VacuumSim::runStep()
{
  runCommonStepTasks();

  initVacuumStep();
  outputVacuumStep();
  runVacuumStep();
}

} /* namespace cosmo */
