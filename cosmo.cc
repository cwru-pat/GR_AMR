#include "cosmo_includes.h"
#include "sims/sim.h"
#include "sims/vacuum.h"
#include "sims/dust.h"
#include "sims/scalar.h"
#include "utils/CartesianCellDoubleQuadraticRefine.h"
#include "utils/CartesianCellDoubleHermiteRefine.h"
#include "utils/CartesianCellDoubleCRSplinesRefine.h"
#include "utils/CartesianCellDoubleCubicRefine.h"
#include "utils/CartesianCellDoubleHermiteCoarsen.h"
#include "utils/CartesianCellDoubleLinearCoarsen.h"
#include "utils/CartesianCellDoubleQuadraticCoarsen.h"
#include "utils/CartesianCellDoubleCRSplinesCoarsen.h"
#include "utils/CartesianCellDoubleCubicCoarsen.h"

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

void add_extra_operators(
  std::shared_ptr<geom::CartesianGridGeometry>& grid_geometry)
{
  grid_geometry->addRefineOperator(
    typeid(pdat::CellVariable<double>).name(),
    std::make_shared<geom::CartesianCellDoubleQuadraticRefine>());
  
  grid_geometry->addRefineOperator(
    typeid(pdat::CellVariable<double>).name(),
    std::make_shared<geom::CartesianCellDoubleHermiteRefine>());

  grid_geometry->addRefineOperator(
    typeid(pdat::CellVariable<double>).name(),
    std::make_shared<geom::CartesianCellDoubleCRSplinesRefine>());

  grid_geometry->addRefineOperator(
    typeid(pdat::CellVariable<double>).name(),
    std::make_shared<geom::CartesianCellDoubleCubicRefine>());

  
  grid_geometry->addCoarsenOperator(
    typeid(pdat::CellVariable<double>).name(),
    std::make_shared<geom::CartesianCellDoubleHermiteCoarsen>());

  grid_geometry->addCoarsenOperator(
    typeid(pdat::CellVariable<double>).name(),
    std::make_shared<geom::CartesianCellDoubleLinearCoarsen>());

  grid_geometry->addCoarsenOperator(
    typeid(pdat::CellVariable<double>).name(),
    std::make_shared<geom::CartesianCellDoubleQuadraticCoarsen>());

  grid_geometry->addCoarsenOperator(
    typeid(pdat::CellVariable<double>).name(),
    std::make_shared<geom::CartesianCellDoubleCRSplinesCoarsen>());

  grid_geometry->addCoarsenOperator(
    typeid(pdat::CellVariable<double>).name(),
    std::make_shared<geom::CartesianCellDoubleCubicCoarsen>());

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

  tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

  
  tbox::pout << "Input file was " << input_filename << std::endl;

  std::string case_name;
  if (argc > 1)
  {
    case_name = argv[1];
  }

  std::shared_ptr<tbox::InputDatabase> input_db(
    new tbox::InputDatabase("input_db"));

  tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

  if (input_db->isDatabase("TimerManager"))
  {
    tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
  }

  std::shared_ptr<tbox::Database> main_db(input_db->getDatabase("Main"));

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

  //set output precision
  int print_precision = main_db->getIntegerWithDefault("print_precision", 9);
  tbox::pout.precision(print_precision);
  tbox::plog.precision(print_precision);

  int num_threads = main_db->getIntegerWithDefault("omp_num_threads", 1);
  // if num_treads == 0, means enable maximum threads
  if(num_threads >= 1)
    omp_set_num_threads(num_threads);

  
  if(main_db->getBoolWithDefault("restart", false))
  {
    std::string restart_name = main_db->getString("restart_basename")+ ".restart";
      
    if(restart_manager->openRestartFile(restart_name,
                                        main_db->getInteger("restart_step"),
                                        main_db->getInteger("restart_nodes")))
    {
      tbox::pout<<"Restarting from file "<<restart_name
                << " with step "<<main_db->getInteger("restart_step")<<"\n";
    }
    else
      tbox::pout<<"Cannot find restart file, will start program from ZERO!\n";
  }
  
  
  std::string simulation_type =
    main_db->getString("simulation_type");

  /*
   * Building hierarchy
   */
  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
    new geom::CartesianGridGeometry(
      dim,
      "CartesianGridGeometry",
      input_db->getDatabase("CartesianGridGeometry")));

  // add extra user defined operators to grid geometry
  add_extra_operators(grid_geometry);
  
  tbox::plog << "Grid Geometry:" << std::endl;
  
  grid_geometry->printClassData(tbox::plog);
  
  std::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
    new hier::PatchHierarchy(
      "Patch Hierarchy",
      grid_geometry,
      input_db->getDatabase("PatchHierarchy")));

  // get file name for VisIt file
  std::string vis_filename =
    main_db->getStringWithDefault("vis_filename", base_name);

  CosmoSim * cosmoSim;

  
  if(simulation_type == "scalar")
  {
    cosmoSim = new ScalarSim(
      patch_hierarchy, dim, input_db, &tbox::plog, simulation_type, vis_filename);
  }
  else if(simulation_type == "vacuum")
  {
    cosmoSim = new VacuumSim(
      patch_hierarchy, dim, input_db, &tbox::plog, simulation_type, vis_filename);
  }
  else if(simulation_type == "dust")
  {
    cosmoSim = new DustSim(
      patch_hierarchy, dim, input_db, &tbox::plog, simulation_type, vis_filename);
  }
  else
  {
    TBOX_ERROR("Invalid simulation type specified.");
  }

  std::shared_ptr<mesh::StandardTagAndInitialize> tag_and_initializer(
    new mesh::StandardTagAndInitialize(
      "CellTaggingMethod",
      cosmoSim,
      input_db->getDatabase("StandardTagAndInitialize")));
  std::shared_ptr<mesh::BergerRigoutsos> box_generator(
    new mesh::BergerRigoutsos(
      dim,
      (input_db->isDatabase("BergerRigoutsos") ?
       input_db->getDatabase("BergerRigoutsos") :
       std::shared_ptr<tbox::Database>())));
  std::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
    new mesh::TreeLoadBalancer(
      dim,
      "load balancer",
      input_db->getDatabase("TreeLoadBalancer")));



  load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());



  std::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
    new mesh::GriddingAlgorithm(
      patch_hierarchy,
      "Gridding Algorithm",
      input_db->getDatabase("GriddingAlgorithm"),
      tag_and_initializer,
      box_generator,
      load_balancer));

  tbox::plog << "Gridding algorithm:" << std::endl;

  //pass the gridding algorithm
  cosmoSim->setGriddingAlgs(gridding_algorithm);

  // Initialize refine and coarsen operators
  cosmoSim->setRefineCoarsenOps(patch_hierarchy);


  // Generate initial conditions
  cosmoSim->setICs(patch_hierarchy);

  // Run simulation
  cosmoSim->run(patch_hierarchy);
  patch_hierarchy->getGridGeometry()->printClassData(tbox::plog);
  
  // print time info
  tbox::TimerManager::getManager()->print(tbox::plog);

  // reset and close stuff
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
