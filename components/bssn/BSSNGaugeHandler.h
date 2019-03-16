/** @file bssn_gauge_fns.h
 * @brief Functions to determine gauge evolution for the BSSN class.
 * Functions are determined via a config setting in the CosmoSim class.
 */

#ifndef COSMO_BSSN_GAUGE_FNS
#define COSMO_BSSN_GAUGE_FNS

#include "../../cosmo_includes.h"
#include "bssn_data.h"

#include <map>

using namespace SAMRAI;

namespace cosmo
{
    
class BSSNGaugeHandler
{
private:
  typedef real_t (BSSNGaugeHandler::*bssn_gauge_func_t)(BSSNData *bd); ///< internal function pointer type

  // Maps to available functions
  std::map<std::string, bssn_gauge_func_t> lapse_gauge_map;
  std::map<std::string, std::map<std::string, bssn_gauge_func_t>> shift_gauge_map;

  // pointers to functions being used
  bssn_gauge_func_t lapse_fn; ///< Lapse evolution function
  bssn_gauge_func_t shift_fn1; ///< Shift evolution function
  bssn_gauge_func_t shift_fn2; ///< Shift evolution function
  bssn_gauge_func_t shift_fn3; ///< Shift evolution function

  // Generic, not evolving gauge
  real_t Static(BSSNData *bd);

  // Harmonic gauge lapse
  real_t HarmonicLapse(BSSNData *bd);
  real_t RelativeHarmonicLapse(BSSNData *bd);

  // 1+log gauge slicing
  real_t gd_c; ///< Tunable gauge parameter
  real_t OnePlusLogLapse(BSSNData *bd);
  real_t RelativeOnePlusLogLapse(BSSNData *bd);
  real_t RelativeAverageOnePlusLogLapse(BSSNData *bd);
  // Untested/experimental lapses
  real_t AnharmonicLapse(BSSNData *bd);
  real_t ConformalSyncLapse(BSSNData *bd);

  // Gamma driver shift function
  real_t GammaDriverShift1(BSSNData *bd);
  real_t GammaDriverShift2(BSSNData *bd);
  real_t GammaDriverShift3(BSSNData *bd);

  // Damped wave gauge
  real_t dw_mu_l; ///< damped wave "mu_l" parameter
  real_t dw_mu_s; ///< damped wave "mu_s" parameter
  real_t dw_p; ///< damped wave "p" parameter
  real_t DampedWaveLapse(BSSNData *bd);
  real_t DampedWaveShift1(BSSNData *bd);
  real_t DampedWaveShift2(BSSNData *bd);
  real_t DampedWaveShift3(BSSNData *bd);

  // AwA Gauge Wave test lapse
  real_t AwAGaugeWaveLapse(BSSNData *bd);

  // AwA Shifted Gauge Wave test gauge
  idx_t AwA_shift_dir; ///< shifted wave direction of prop. (\in {1,2,3})
  real_t AwAShiftedWaveLapse(BSSNData *bd);
  real_t AwAShiftedWaveShift1(BSSNData *bd);
  real_t AwAShiftedWaveShift2(BSSNData *bd);
  real_t AwAShiftedWaveShift3(BSSNData *bd);

  // Map of strings to functions
  void _initGaugeMaps()
  {
    // Lapse functions
    lapse_gauge_map["Static"] = &BSSNGaugeHandler::Static;
    lapse_gauge_map["Harmonic"] = &BSSNGaugeHandler::HarmonicLapse;
    lapse_gauge_map["RelativeHarmonic"] = &BSSNGaugeHandler::RelativeHarmonicLapse;
    lapse_gauge_map["Anharmonic"] = &BSSNGaugeHandler::AnharmonicLapse;
    lapse_gauge_map["OnePlusLog"] = &BSSNGaugeHandler::OnePlusLogLapse;
    lapse_gauge_map["RelativeOnePlusLog"] = &BSSNGaugeHandler::RelativeOnePlusLogLapse;
    lapse_gauge_map["RelativeAverageOnePlusLog"] = &BSSNGaugeHandler::RelativeAverageOnePlusLogLapse;
    lapse_gauge_map["DampedWave"] = &BSSNGaugeHandler::DampedWaveLapse;
    lapse_gauge_map["ConformalSync"] = &BSSNGaugeHandler::ConformalSyncLapse;
    lapse_gauge_map["AwAGaugeWave"] = &BSSNGaugeHandler::AwAGaugeWaveLapse;
    lapse_gauge_map["AwAShiftedWave"] = &BSSNGaugeHandler::AwAShiftedWaveLapse;

    // Shift functions
    // Static gauge
    shift_gauge_map["Static"]["1"] = &BSSNGaugeHandler::Static;
    shift_gauge_map["Static"]["2"] = &BSSNGaugeHandler::Static;
    shift_gauge_map["Static"]["3"] = &BSSNGaugeHandler::Static;
    // gamma driver
    shift_gauge_map["GammaDriver"]["1"] = &BSSNGaugeHandler::GammaDriverShift1;
    shift_gauge_map["GammaDriver"]["2"] = &BSSNGaugeHandler::GammaDriverShift2;
    shift_gauge_map["GammaDriver"]["3"] = &BSSNGaugeHandler::GammaDriverShift3;
    // Damped wave
    shift_gauge_map["DampedWave"]["1"] = &BSSNGaugeHandler::DampedWaveShift1;
    shift_gauge_map["DampedWave"]["2"] = &BSSNGaugeHandler::DampedWaveShift2;
    shift_gauge_map["DampedWave"]["3"] = &BSSNGaugeHandler::DampedWaveShift3;
    // AwA shifted wave test
    shift_gauge_map["AwAShiftedWave"]["1"] = &BSSNGaugeHandler::AwAShiftedWaveShift1;
    shift_gauge_map["AwAShiftedWave"]["2"] = &BSSNGaugeHandler::AwAShiftedWaveShift2;
    shift_gauge_map["AwAShiftedWave"]["3"] = &BSSNGaugeHandler::AwAShiftedWaveShift3;
  }

  void _initDefaultParameters(std::shared_ptr<tbox::Database> database)
  {
    AwA_shift_dir = database->getIntegerWithDefault("AwA_shift_dir", 1);
   
    dw_mu_l = database->getDoubleWithDefault("dw_mu_l", 0.0);
    dw_mu_s = database->getDoubleWithDefault("dw_mu_s", 0.0);
    dw_p = database->getDoubleWithDefault("dw_p", 0.0);
    gd_c = database->getDoubleWithDefault("gd_c", 1.0);
  }

public:

  /**
   * @brief Initialize with static, non-evolving gauge
   */


  /**
   * @brief Initialize with gauge determined by config file (default to a "static", non-evolving gauge)
   */
  BSSNGaugeHandler(std::shared_ptr<tbox::Database> database)
  {
    _initGaugeMaps();
    _initDefaultParameters(database);
    setLapseFn(database->getStringWithDefault("lapse", "Static"));
    setShiftFn(database->getStringWithDefault("Shift", "Static"));
  }

  /**
   * @brief Set the lapse function
   */
  void setLapseFn(std::string name)
  {
    if ( lapse_gauge_map.find(name) == lapse_gauge_map.end() )
    {
      TBOX_ERROR("Error: Lapse gauge not found: `" << name << "`!\n");
    }

    tbox::plog<<"Setting lapse function with "<<name<<"\n";
    lapse_fn = lapse_gauge_map[name];
  }

  /**
   * @brief Set the shift function
   */
  void setShiftFn(std::string name)
  {

    tbox::plog<<"Setting shift function with "<<name<<"\n";
    shift_fn1 = shift_gauge_map[name]["1"];
    shift_fn2 = shift_gauge_map[name]["2"];
    shift_fn3 = shift_gauge_map[name]["3"];
  }

  /**
   * @brief Lapse evolution function for BSSN class to call
   */
  real_t ev_lapse(BSSNData *bd)
  {
    return (*this.*lapse_fn)(bd);
  }

  /**
   * @brief Shift in x-dir evolution function for BSSN class to call
   */
  real_t ev_shift1(BSSNData *bd)
  {
    return (*this.*shift_fn1)(bd);
  }

  /**
   * @brief Shift in y-dir evolution function for BSSN class to call
   */
  real_t ev_shift2(BSSNData *bd)
  {
    return (*this.*shift_fn2)(bd);
  }

  /**
   * @brief Shift in z-dir evolution function for BSSN class to call
   */
  real_t ev_shift3(BSSNData *bd)
  {
    return (*this.*shift_fn3)(bd);
  }

};

}

#endif
