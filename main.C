#include "SAMRAI/SAMRAI_config.h"

#include IOMANIP_HEADER_FILE
#include <fstream>
// using namespace std;

#include "WaveSolv.h"
#include "SAMRAI/source/test/testlib/get-input-filename.h"

/*
 * Headers for basic SAMRAI objects used in this code.
 */
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/TimeInterpolateOperator.h"
#include "SAMRAI/xfer/PatchLevelFillPattern.h"
#include "SAMRAI/xfer/PatchLevelBorderFillPattern.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"
#include "SAMRAI/xfer/PatchLevelInteriorFillPattern.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/geom/CartesianCellDoubleLinearRefine.h"

/*
 * Headers for major algorithm/data structure objects from SAMRAI
 */
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/solv/FACPreconditioner.h"

#include "boost/shared_ptr.hpp"

#include <vector>
#include <string>

using namespace SAMRAI;

int main(int argc, char* argv[])
{
  std::string input_filename;

  /*
   * Initialize MPI, process argv, and initialize SAMRAI
   */
  tbox::SAMRAI_MPI::init(&argc, &argv);
  if (get_input_filename(&argc, argv, input_filename) == 1) {
    tbox::pout << "Usage: " << argv[0] << " <input file>." << std::endl;
    tbox::SAMRAI_MPI::finalize();
    return 0;
  }
  tbox::SAMRAIManager::initialize();
  tbox::SAMRAIManager::startup();


  tbox::pout << "Input file is " << input_filename << std::endl;

  std::string case_name;
  if (argc > 1)
  {
    case_name = argv[1];
  }

  boost::shared_ptr<tbox::InputDatabase> input_db(
    new tbox::InputDatabase("input_db"));
  tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

  if (input_db->isDatabase("TimerManager"))
  {
    tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
  }


  boost::shared_ptr<tbox::Database> main_db(input_db->getDatabase("Main"));

  const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

  tbox::plog << "Main database:" << std::endl;
  main_db->printClassData(tbox::plog);

  std::string base_name =
    main_db->getStringWithDefault("base_name", "noname");

  
  if (!case_name.empty()) {
    base_name = base_name + '-' + case_name;
  }
  base_name = base_name + '-'
    + tbox::Utilities::intToString(tbox::SAMRAI_MPI::getSAMRAIWorld().getSize(), 5);

  /*
   * Log file info.
   */
  std::string log_filename =
    main_db->getStringWithDefault("log_filename", base_name + ".log");
  bool log_all =
    main_db->getBoolWithDefault("log_all", false);
  if (log_all)
    tbox::PIO::logAllNodes(log_filename);
  else
    tbox::PIO::logOnlyNodeZero(log_filename);


  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
    new geom::CartesianGridGeometry(
      dim,
      "CartesianGridGeometry",
      input_db->getDatabase("CartesianGridGeometry")));

  tbox::plog << "Grid Geometry:" << std::endl;
  grid_geometry->printClassData(tbox::plog);
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
    new hier::PatchHierarchy(
      "Patch Hierarchy",
      grid_geometry,
      input_db->getDatabase("PatchHierarchy")));
  
  WaveSolv solver(dim,
                  *(input_db->getDatabase("WaveSolv")),
                  &tbox::plog);



  
  boost::shared_ptr<mesh::StandardTagAndInitialize> tag_and_initializer(
    new mesh::StandardTagAndInitialize(
      "CellTaggingMethod",
      &solver,
      input_db->getDatabase("StandardTagAndInitialize")));
  boost::shared_ptr<mesh::BergerRigoutsos> box_generator(
    new mesh::BergerRigoutsos(
      dim,
      (input_db->isDatabase("BergerRigoutsos") ?
       input_db->getDatabase("BergerRigoutsos") :
       boost::shared_ptr<tbox::Database>())));
  boost::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
    new mesh::TreeLoadBalancer(
      dim,
      "load balancer",
      input_db->getDatabase("TreeLoadBalancer")));

  load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

  boost::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
    new mesh::GriddingAlgorithm(
      patch_hierarchy,
      "Gridding Algorithm",
      input_db->getDatabase("GriddingAlgorithm"),
      tag_and_initializer,
      box_generator,
      load_balancer));

  tbox::plog << "Gridding algorithm:" << std::endl;
  gridding_algorithm->printClassData(tbox::plog);
  /*
   * Make the coarse patch level.
   */

  gridding_algorithm->makeCoarsestLevel(0.0);

  int ln;

  
  double dt = grid_geometry->getDx()[0] / 10.0, cur_t = 0.0, fin_t;

  int step_cnt = 0;
  
  fin_t = main_db->getDouble("final_time");

  //  throw(-1);
  
  while(cur_t < fin_t)
  {    
    
    /*for(int ln = 0; ln < 1; ln ++)
    {
      boost::shared_ptr<hier::PatchLevel> level(patch_hierarchy->getPatchLevel(ln));   
         for (hier::PatchLevel::iterator pi(level->begin());
        pi != level->end(); ++pi) {
      hier::Patch& patch = **pi;
      boost::shared_ptr<pdat::CellData<double> > current_solution(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(solver.d_phi_previous)));
      const hier::Box& ghost_box = current_solution->getGhostBox();
      current_solution->print(ghost_box, 0, std::cout, 4);
      }
      }
      std::cout<<"\n\n\n\n\n";*/
    if(step_cnt % 2 == 0)
    {
      std::vector<int> tag_buffer(patch_hierarchy->getMaxNumberOfLevels());
      for (ln = 0; ln < static_cast<int>(tag_buffer.size()); ++ln)
      {
        tag_buffer[ln] = 1;
      }
      //if(step_cnt == 0)
      gridding_algorithm->regridAllFinerLevels(
        0,
        tag_buffer,
        0,
        0.0);
      
      tbox::plog << "Newly adapted hierarchy\n";
      patch_hierarchy->recursivePrint(tbox::plog, "    ", 1);
    }


    
    solver.advanceHierarchy(patch_hierarchy, cur_t, cur_t + dt);
    cur_t += dt;
    step_cnt++;
        std::cout<<"Max error is "<<solver._maxError(patch_hierarchy)<<"\n";
    /*for(int ln = 0; ln < 2; ln ++)
    {
      boost::shared_ptr<hier::PatchLevel> level(patch_hierarchy->getPatchLevel(ln));   
         for (hier::PatchLevel::iterator pi(level->begin());
        pi != level->end(); ++pi) {
      hier::Patch& patch = **pi;
      boost::shared_ptr<pdat::CellData<double> > current_solution(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(solver.d_phi_previous)));
      const hier::Box& ghost_box = current_solution->getGhostBox();
      current_solution->print(ghost_box, 0, std::cout, 4);
      }
      }*/

    
    //      throw(-1);
  }
  tbox::TimerManager::getManager()->print(tbox::plog);
  tbox::SAMRAIManager::shutdown();
  tbox::SAMRAIManager::finalize();
  tbox::SAMRAI_MPI::finalize();

}
