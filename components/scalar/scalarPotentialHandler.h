#ifndef COSMO_SCALAR_POTENTIAL_FNS
#define COSMO_SCALAR_POTENTIAL_FNS

#include "../../cosmo_includes.h"
#include "scalar_data.h"
#include "../components/bssn/bssn_data.h"
#include <map>

using namespace SAMRAI;

namespace cosmo
{
  
class scalarPotentialHandler
{
private:
  typedef real_t (scalarPotentialHandler::*scalar_potential_func_t)(
    BSSNData *bd, ScalarData *sd); ///< internal function pointer type

  // Maps to available functions
  std::map<std::string, scalar_potential_func_t> scalar_potential_map;
  std::map<std::string, scalar_potential_func_t> scalar_der_potential_map;
  
  // pointers to functions being used
  scalar_potential_func_t potential_fn; ///< Lapse evolution function
  scalar_potential_func_t der_potential_fn; ///< Lapse evolution function

  // constant potential
  real_t constant(BSSNData *bd, ScalarData *sd);
  // constant potential
  real_t der_constant(BSSNData *bd, ScalarData *sd);

  // 
  real_t quadratic(BSSNData *bd, ScalarData *sd);
  // 
  real_t der_quadratic(BSSNData *bd, ScalarData *sd);

  real_t Lambda, q_coef; 

  // Map of strings to functions
  void _initGaugeMaps()
  {
    // Lapse functions
    scalar_potential_map["Constant"] = &scalarPotentialHandler::constant;
    scalar_potential_map["Quadratic"] = &scalarPotentialHandler::quadratic;
    scalar_der_potential_map["Constant"] = &scalarPotentialHandler::der_constant;
    scalar_der_potential_map["Quadratic"] = &scalarPotentialHandler::der_quadratic;

  }

  void _initDefaultParameters(boost::shared_ptr<tbox::Database> database)
  {
    Lambda = database->getDoubleWithDefault("Lambda", 0.0);
   
    q_coef = database->getDoubleWithDefault("q_coef", 0.0);

  }

public:

  /**
   * @brief Initialize with static, non-evolving gauge
   */


  /**
   * @brief Initialize with gauge determined by config file (default to a "static", non-evolving gauge)
   */
  scalarPotentialHandler(boost::shared_ptr<tbox::Database> database)
  {
    _initGaugeMaps();
    _initDefaultParameters(database);
    setPotential(database->getStringWithDefault("potential_type", "Constant"));
  }

  /**
   * @brief Set the lapse function
   */
  void setPotential(std::string name)
  {
    if (scalar_potential_map.find(name) == scalar_potential_map.end())
    {
      TBOX_ERROR("Error: Lapse gauge not found: `" << name << "`!\n");
    }

    tbox::plog<<"Setting lapse function with "<<name<<"\n";
    potential_fn = scalar_potential_map[name];
    der_potential_fn = scalar_der_potential_map[name];
  }


  /**
   * @brief Lapse evolution function for BSSN class to call
   */
  real_t ev_potential(BSSNData *bd, ScalarData *sd)
  {
    return (*this.*potential_fn)(bd, sd);
  }
  real_t ev_der_potential(BSSNData *bd, ScalarData *sd)
  {
    return (*this.*der_potential_fn)(bd, sd);
  }

};

}

#endif
