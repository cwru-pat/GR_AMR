/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   AdaptivePoisson class implementation
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/pdat/ArrayDataAccess.h"
#include "SAMRAI/solv/CellPoissonFACOps.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/geom/CartesianCellDoubleLinearRefine.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/geom/CartesianCellDoubleLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianSideDoubleWeightedAverage.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

#include "boost/shared_ptr.hpp"
#include "WaveSolv.h"
#include <sstream>
#include <iomanip>
#include <cstring>
#include <stdlib.h>

#include <cmath>

using namespace SAMRAI;

#define PI 3.14159265359

#define PW2(x) ((x) * (x))

boost::shared_ptr<tbox::Timer> WaveSolv::t_advance_hier;

WaveSolv::WaveSolv(
  const tbox::Dimension& dim,
  tbox::Database& database,
  std::ostream* log_stream = 0):
  d_lstream(log_stream),
  d_database(&database),
  d_barrier_and_time(false),
  d_dim(dim),
  d_omega(0),
  d_context_current(new hier::VariableContext("CURRENT")),
  d_context_scratch(new hier::VariableContext("SCRATCH")),
  d_context_previous(new hier::VariableContext("PREVIOUS")),
  d_phi(new pdat::CellVariable<double>(d_dim, "var:phi", 1)),
  d_pi(new pdat::CellVariable<double>(d_dim, "var:pi", 1)),
  d_wt(new pdat::CellVariable<double>(d_dim, "var:wt", 1)),
  nx(database.getDouble("nx")),
  ny(database.getDouble("ny")),
  nz(database.getDouble("nz")),
  d_adaption_threshold(2)
{
  
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  d_phi_current =
    variable_db->registerVariableAndContext(
      d_phi,
      d_context_current,
      hier::IntVector(d_dim, 1));
  d_phi_scratch =
    variable_db->registerVariableAndContext(
      d_phi,
      d_context_scratch,
      hier::IntVector(d_dim, 1));
  d_phi_previous =
    variable_db->registerVariableAndContext(
      d_phi,
      d_context_previous,
      hier::IntVector(d_dim, 1));

  d_pi_current =
    variable_db->registerVariableAndContext(
      d_pi,
      d_context_current,
      hier::IntVector(d_dim, 1));
  d_pi_scratch =
    variable_db->registerVariableAndContext(
      d_pi,
      d_context_scratch,
      hier::IntVector(d_dim, 1));
  d_pi_previous =
    variable_db->registerVariableAndContext(
      d_pi,
      d_context_previous,
      hier::IntVector(d_dim, 1));
  d_weight =
    variable_db->registerVariableAndContext(
      d_wt,
      d_context_current,
      hier::IntVector(d_dim, 0));
}


double WaveSolv::_laplacian(
  boost::shared_ptr<pdat::CellData<double>>& var,
    int i, int j, int k, const double dx[])
{
  MDA_Access<double, 3, MDA_OrderColMajor<3> > va =
    pdat::ArrayDataAccess::access<3, double>(var->getArrayData());

  return (va(i+1,j,k) + va(i-1,j,k) - 2.0 * va(i,j,k)) / dx[0] / dx[0]
    + (va(i,j+1,k) + va(i,j-1,k) - 2.0 * va(i,j,k)) / dx[1] / dx[1]
    + (va(i,j,k+1) + va(i,j,k-1) - 2.0 * va(i,j,k)) / dx[2] / dx[2];
}

double WaveSolv::_der_norm(
  boost::shared_ptr<pdat::CellData<double>>& var,
    int i, int j, int k, const double dx[])
{
  MDA_Access<double, 3, MDA_OrderColMajor<3> > va =
    pdat::ArrayDataAccess::access<3, double>(var->getArrayData());

  return sqrt(PW2(va(i+1,j,k) + va(i-1,j,k))
             + PW2(va(i,j+1,k) + va(i,j-1,k))
              + PW2(va(i,j,k+1) + va(i,j,k-1)));
}


double WaveSolv::_maxError(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  double max_err = 0;
  
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
      
      boost::shared_ptr<pdat::CellData<double> > phi(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(d_phi_previous)));

      boost::shared_ptr<pdat::CellData<double> > weight(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(d_weight)));
      
      MDA_AccessConst<double, 3, MDA_OrderColMajor<3> > phi_array =
        pdat::ArrayDataAccess::access<3, double>(
          phi->getArrayData());

       MDA_AccessConst<double, 3, MDA_OrderColMajor<3> > weight_array =
        pdat::ArrayDataAccess::access<3, double>(
          weight->getArrayData());

       const double * domain_lower = &grid_geometry.getXLower()[0];
       const double * domain_upper = &grid_geometry.getXUpper()[0];
       
      
      
       double L[3];

       for(int i = 0 ; i < 3; i++)
         L[i] = domain_upper[i] - domain_lower[i];

       
       const int * lower = &box.lower()[0];
       const int * upper = &box.upper()[0];
      

       double t = patch->getPatchData(d_phi_previous)->getTime();
       //       std::cout<<t<<"\n";
       const double *dx = &patch_geom->getDx()[0];
       for(int k = lower[2]; k <= upper[2]; k++)
       {
         for(int j = lower[1]; j <= upper[1]; j++)
         {
           for(int i = lower[0]; i <= upper[0]; i++)
           {
             if(weight_array(i,j,k) > 0)
             {
               max_err = tbox::MathUtilities<double>::Max(
                 max_err, tbox::MathUtilities<double>::Abs(
                   sin(2.0 * PI * (double) i * dx[0] * nx / L[0] +
                       2.0 * PI * (double) j * dx[1] * ny / L[1] +
                       2.0 * PI * (double) k * dx[2] * nz / L[2] +
                       d_omega * t) -
                   phi_array(i,j,k)));
               /*                   if(t>=0.09999)
                   {
                     std::cout<<tbox::MathUtilities<double>::Abs(
                   sin(2.0 * PI * (double) i * dx[0] * nx / L[0] +
                       2.0 * PI * (double) j * dx[1] * ny / L[1] +
                       2.0 * PI * (double) k * dx[2] * nz / L[2] +
                       d_omega * t) - phi_array(i,j,k))<<" "<<max_err<<"\n";
                       };*/
             }

           }
         }
       }

    }
  }
  return max_err;
}

void WaveSolv::initCoarsest(
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
    //double nx = d_database->getDouble("nx");
    //double ny = d_database->getDouble("ny");
    //double nz = d_database->getDouble("nz");

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

void WaveSolv::computeVectorWeights(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  TBOX_ASSERT(hierarchy);
  TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *hierarchy);

  int weight_id = d_weight;
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

void WaveSolv::initializeLevelData(
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
     level->allocatePatchData(d_phi_current);
     level->allocatePatchData(d_phi_previous);
     level->allocatePatchData(d_pi_current);
     level->allocatePatchData(d_pi_previous);
     level->allocatePatchData(d_weight);
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
       grid_geometry.
       lookupRefineOperator(d_phi, "LINEAR_REFINE");
     
          TBOX_ASSERT(accurate_refine_op);
     

     refiner.registerRefine(d_phi_previous,
         d_phi_previous,
         d_phi_previous,
         accurate_refine_op);

     refiner.registerRefine(d_pi_previous,
         d_pi_previous,
         d_pi_previous,
         accurate_refine_op);
      
      boost::shared_ptr<xfer::RefineSchedule> refine_schedule;


      
      if (ln > 0) {
         /*
          * Include coarser levels in setting data
          */
         refine_schedule =
            refiner.createSchedule(level,
               old_level,
               ln - 1,
               hierarchy,
               NULL);
      } else {
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

      math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);
      if (refine_schedule) {
         refine_schedule->fillData(0.0);
         // It is null if this is the bottom level.
      } else {
        //hcellmath.setToScalar(d_phi_current, 0.0, false);
        //hcellmath.setToScalar(d_pi_current, 0.0, false);
        std::cout<<"Can not get refine schedul, check your code!\n";
        throw(-1);
      }

      if (0) {
         // begin debug code
         math::HierarchyCellDataOpsReal<double> hcellmath_debug(hierarchy);
         hcellmath_debug.printData(d_phi_current, tbox::pout, false);
         // end debug code
      }
   }

   /* Set vector weight. */
   computeVectorWeights(hierarchy);
}

void WaveSolv::resetHierarchyConfiguration(
   /*! New hierarchy */ const boost::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
   /*! Coarsest level */ int coarsest_level,
   /*! Finest level */ int finest_level)
{
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);

   d_hierarchy = new_hierarchy;
   /*
    * Recompute or reset internal data tied to the hierarchy,
    * if any.  None at this time.
    */
   /*
    * Log the new hierarchy.
    */
   if (d_lstream) {
      *d_lstream
      << "AdaptivePoisson::resetHierarchyConfiguration\n";
      d_hierarchy->recursivePrint(*d_lstream, "    ", 2);
   }
}



void WaveSolv::applyGradientDetector(
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

   if (d_lstream) {
      *d_lstream
      << "WaveSolv("  << ")::applyGradientDetector"
      << std::endl;
   }
   hier::PatchHierarchy& hierarchy = *hierarchy_;
   boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
     BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
       hierarchy.getGridGeometry()));
   double max_lap = 0;
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
      
      boost::shared_ptr<pdat::CellData<double>> phi_data(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(patch.getPatchData(d_phi_previous)));

      boost::shared_ptr<pdat::CellData<double> > weight_(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
          patch.getPatchData(d_weight)));

      if (!phi_data) {
         TBOX_ERROR("Data index " << d_phi_current
                                  << " does not exist for patch.\n");
      }
      pdat::CellData<int>& tag_cell_data = *tag_cell_data_;
      pdat::CellData<double> & weight = *weight_;
          
      tag_cell_data.fill(0);
      
      hier::Box::iterator iend(patch.getBox().end());

      for (hier::Box::iterator i(patch.getBox().begin()); i != iend; ++i)
      {

         const pdat::CellIndex cell_index(*i);
         //if(ln == 0 &&  hierarchy.getNumberOfLevels() > 1) std::cout<<weight(cell_index)<<"\n";
         /*if( tbox::MathUtilities<double>::Abs(
              _laplacian(
              phi_data,
              cell_index(0),
              cell_index(1),
              cell_index(2),
              &(patch_geom->getDx())[0])) > d_adaption_threshold)*/
         if(_der_norm(
              phi_data,
              cell_index(0),
              cell_index(1),
              cell_index(2),
              &(patch_geom->getDx())[0]) > d_adaption_threshold)
         {
           /*max_lap = tbox::MathUtilities<double>::Max(
             max_lap, tbox::MathUtilities<double>::Abs(
             _laplacian(
               phi_data,
               cell_index(0),
               cell_index(1),
               cell_index(2),
               &patch_geom->getDx()[0])));*/
           max_lap = _der_norm(
              phi_data,
              cell_index(0),
              cell_index(1),
              cell_index(2),
              &(patch_geom->getDx())[0]);
           tag_cell_data(cell_index) = 1;
           ++ntag;
           
         }
         
      }

   }
   tbox::plog << "Adaption threshold is " << d_adaption_threshold << "\n";
   tbox::plog << "Number of cells tagged on level " << ln << " is "
              << ntag << "/" << ntotal << "\n";
   tbox::plog << "Max norm is " << max_lap << "\n";
}

void WaveSolv::advanceLevel(
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


  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;

    const hier::Box& box = patch->getBox();

    boost::shared_ptr<pdat::CellData<double> > pi_current(
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
        patch->getPatchData(d_pi_current)));
    boost::shared_ptr<pdat::CellData<double> > phi_current(
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
        patch->getPatchData(d_phi_current)));
    boost::shared_ptr<pdat::CellData<double> > pi_previous(
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
        patch->getPatchData(d_pi_previous)));
    boost::shared_ptr<pdat::CellData<double> > phi_previous(
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
        patch->getPatchData(d_phi_previous)));


    MDA_Access<double, 3, MDA_OrderColMajor<3> > pi_current_array =
      pdat::ArrayDataAccess::access<3, double>(
        pi_current->getArrayData());

    MDA_Access<double, 3, MDA_OrderColMajor<3> > phi_current_array =
      pdat::ArrayDataAccess::access<3, double>(
        phi_current->getArrayData());
    MDA_Access<double, 3, MDA_OrderColMajor<3> > pi_previous_array =
      pdat::ArrayDataAccess::access<3, double>(
        pi_previous->getArrayData());

    MDA_Access<double, 3, MDA_OrderColMajor<3> > phi_previous_array =
      pdat::ArrayDataAccess::access<3, double>(
        phi_previous->getArrayData());

    patch->getPatchData(d_phi_current)->setTime(to_t);
    patch->getPatchData(d_pi_current)->setTime(to_t);


    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry()));

    //
    
    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {

          phi_current_array(i,j,k) =
            phi_previous_array(i,j,k) + pi_previous_array(i,j,k) * dt;
          
          pi_current_array(i,j,k) =
            pi_previous_array(i,j,k) +
            _laplacian(phi_previous, i,j,k, &patch_geom->getDx()[0]) * dt;
                   if(from_t <= 0.05 && from_t >= 0.04 )
          {
            //std::cout<<phi_current_array(i,j,k)<<" "<<_laplacian(phi_previous, i,j,k, &patch_geom->getDx()[0])<<" "<<
              //pi_current_array(i,j,k)<<"\n";
          }
        }
      }
    }

  }

  xfer::RefineAlgorithm refiner;
  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  boost::shared_ptr<hier::RefineOperator> space_refine_op =
       grid_geometry.
       lookupRefineOperator(d_phi, "CONSERVATIVE_LINEAR_REFINE");
  boost::shared_ptr<hier::TimeInterpolateOperator> time_refine_op =
    grid_geometry.
    lookupTimeInterpolateOperator(d_phi, "STD_LINEAR_TIME_INTERPOLATE");
  boost::shared_ptr<xfer::RefineSchedule> refine_schedule;

  
  
  if(ln ==0)
  {
    refiner.registerRefine(d_phi_current,
                           d_phi_current,
                           d_phi_current,
                           space_refine_op);
    refiner.registerRefine(d_pi_current,
                           d_pi_current,
                           d_pi_current,
                           space_refine_op);

    refine_schedule = refiner.createSchedule(level, NULL);

    refine_schedule->fillData(to_t);
  }
  else
  {
    refiner.registerRefine(d_phi_current,
                           d_phi_current,
                           d_phi_previous,
                           d_phi_current,
                           d_phi_scratch,
                           space_refine_op,
                           time_refine_op);
    refiner.registerRefine(d_phi_current,
                           d_phi_current,
                           d_phi_previous,
                           d_phi_current,
                           d_phi_scratch,
                           space_refine_op,
                           time_refine_op);

    refine_schedule = refiner.createSchedule(level, level, ln-1, hierarchy, NULL, true, NULL);
    
    refine_schedule->fillData(to_t);
  }

  advanceLevel(hierarchy, ln+1, from_t, from_t + (to_t - from_t)/2.0);

  advanceLevel(hierarchy, ln+1, from_t + (to_t - from_t)/2.0, to_t);

  if(ln < hierarchy->getNumberOfLevels() -1 )
  {
    xfer::CoarsenAlgorithm coarsener(d_dim);

    boost::shared_ptr<hier::CoarsenOperator> space_coarsen_op =
      grid_geometry.
      lookupCoarsenOperator(d_phi, "CONSERVATIVE_COARSEN");
  

    boost::shared_ptr<xfer::CoarsenSchedule> coarsen_schedule;

    coarsener.registerCoarsen(d_phi_current,
                              d_phi_current,
                              space_coarsen_op,
                              NULL);
  
    coarsener.registerCoarsen(d_pi_current,
                              d_pi_current,
                              space_coarsen_op,
                              NULL);
  
    coarsen_schedule = coarsener.createSchedule(level, hierarchy->getPatchLevel(ln+1));
    coarsen_schedule->coarsenData();

    /*    xfer::RefineAlgorithm refiner_after;
    refiner_after.registerRefine(d_phi_current,
                           d_phi_current,
                           d_phi_current,
                           space_refine_op);
    refiner_after.registerRefine(d_pi_current,
                           d_pi_current,
                           d_pi_current,
                           space_refine_op);

    refine_schedule = refiner_after.createSchedule(level, NULL);

    refine_schedule->fillData(to_t);*/

  }

  math::PatchCellDataOpsReal<double> pcellmath;

  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;

    pcellmath.swapData(patch, d_phi_current, d_phi_previous);
    pcellmath.swapData(patch, d_pi_current, d_pi_previous);
  }
  
}

void WaveSolv::advanceHierarchy(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  double from_t, double to_t)
{
  //  t_advance_hier->start();
  if(d_barrier_and_time)
    t_advance_hier->barrierAndStart();

  advanceLevel(hierarchy,
               0,
               from_t,
               to_t);
  if (d_barrier_and_time)
  {
    t_advance_hier->stop();
  }

}
