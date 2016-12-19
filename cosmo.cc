#include "cosmo_includes.h"
#include "sims/sim.h"
#include "sims/vacuum.h"

using namespace SAMRAI;
using namespace cosmo;

int get_input_filename(
   int* argc,
   char* argv[],
   std::string& input_filename) {

   int rval = 0;
   std::string argv0(argv[0]);
   if (*argc > 1) {
      // Input file is the first argument.  Shift other arguments down.
      input_filename = argv[1];
      --(*argc);
      int i;
      for (i = 1; i < (*argc); ++i) {
         argv[i] = argv[i + 1];
      }
   } else if (*argc == 1 && argv0.rfind("check-") < argv0.size()) {
      /*
       * No argument but input file is implicit in program name
       * which has the form check-<input file without .input>.
       */
      input_filename = argv0.substr(argv0.rfind("check-") + 6) + ".input";
   } else if (*argc == 1) {
      // No argument and not invoked as "check-blah".
      rval = 1;
   }
   return rval;
}



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


  tbox::pout << "Input file was " << input_filename << std::endl;

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

  tbox::pout.precision(PRINT_PRECISION);
  tbox::plog.precision(PRINT_PRECISION);

  
  std::string simulation_type =
    main_db->getString("simulation_type");

  /*
   * Building hierarchy
   */
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
 
 

  

  std::string vis_filename =
    main_db->getStringWithDefault("vis_filename", base_name);


  CosmoSim * cosmoSim;

  
  if(simulation_type == "scalar")
  {
    //cosmoSim = new ScalarSim();
    TBOX_ERROR("Have not finished scalar part yet.");
  }
  else if(simulation_type == "vacuum")
  {
    cosmoSim = new VacuumSim(
      dim, input_db, &tbox::plog, simulation_type, vis_filename);
  }
  else
  {
    TBOX_ERROR("Invalid simulation type specified.");
  }

  boost::shared_ptr<mesh::StandardTagAndInitialize> tag_and_initializer(
    new mesh::StandardTagAndInitialize(
      "CellTaggingMethod",
      cosmoSim,
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



  //deliver the gridding algorithm
  cosmoSim->setGriddingAlgs(gridding_algorithm);

  // Initialize simulation 

  cosmoSim->setRefineCoarsenOps(patch_hierarchy);
  // Generate initial conditions
  cosmoSim->setICs();



  cosmoSim->setRefineCoarsenOps(patch_hierarchy);
  // Run simulation
  cosmoSim->run(patch_hierarchy);

  
  
  tbox::TimerManager::getManager()->print(tbox::plog);
  
  patch_hierarchy.reset();
  load_balancer.reset();
  gridding_algorithm.reset();
  grid_geometry.reset();
  box_generator.reset();
  tbox::SAMRAIManager::shutdown();
  tbox::SAMRAIManager::finalize();
  tbox::SAMRAI_MPI::finalize();

  return 0;
}
