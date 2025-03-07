/**
 * @file   geometry_sbnd_test.cxx
 * @brief  Unit test for geometry on SBND detector
 * @date   May 11th, 2015
 * @author petrillo@fnal.gov
 * 
 * Usage:
 *   `geometry_sbnd_test  [ConfigurationFile [GeometryTestParameterSet]]`
 * 
 * By default, GeometryTestParameterSet is set to `"physics.analyzers.geotest"`.
 * 
 */

// SBND libraries
#include "sbndcode/Geometry/GeoObjectSorterSBND.h"
#include "sbndcode/Geometry/WireReadoutSorterSBND.h"

// LArSoft libraries
#include "test/Geometry/geometry_unit_test_sbnd.h"
#include "larcorealg/test/Geometry/GeometryTestAlg.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larcorealg/Geometry/AuxDetReadoutGeom.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutStandardGeom.h"
#include "larcorealg/Geometry/StandaloneGeometrySetup.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ standard libraries
#include <memory>


//------------------------------------------------------------------------------
//---  The test environment
//---


// we define here all the configuration that is needed;
// we use an existing class provided for this purpose, since our test
// environment allows us to tailor it at run time.
using SBNDGeometryConfiguration
  = sbnd::testing::SBNDGeometryEnvironmentConfiguration;

/*
 * GeometryTesterFixture, configured with the object above, is used in a
 * non-Boost-unit-test context.
 * It provides:
 * - `geo::GeometryCore const* Geometry()`
 * - `geo::GeometryCore const* GlobalGeometry()` (static member)
 */
using SBNDGeometryTestEnvironment
  = testing::GeometryTesterEnvironment<SBNDGeometryConfiguration, geo::GeoObjectSorterSBND>;


//------------------------------------------------------------------------------
//---  The tests
//---

/** ****************************************************************************
 * @brief Runs the test
 * @param argc number of arguments in argv
 * @param argv arguments to the function
 * @return number of detected errors (0 on success)
 * @throw cet::exception most of error situations throw
 * 
 * The arguments in argv are:
 * 0. name of the executable ("Geometry_test")
 * 1. path to the FHiCL configuration file
 * 2. FHiCL path to the configuration of the geometry test
 *    (default: physics.analysers.geotest)
 * 3. FHiCL path to the configuration of the geometry
 *    (default: services.Geometry)
 * 
 */
//------------------------------------------------------------------------------
int main(int argc, char const** argv) {
  
  SBNDGeometryConfiguration config("geometry_test_SBND");
  config.SetMainTesterParameterSetName("geotest");
  
  //
  // parameter parsing
  //
  int iParam = 0;
  
  // first argument: configuration file (mandatory)
  if (++iParam < argc) config.SetConfigurationPath(argv[iParam]);
  
  // second argument: path of the parameter set for geometry test configuration
  // (optional; default: "physics.analysers.geotest");
  // if no path is provided, we have a empty default configuration;
  // if path is provided, we don't have any default configuration
  // and if the configuration is missing there will be an error
  if (++iParam < argc) config.SetMainTesterParameterSetPath(argv[iParam]);
  else                 config.AddDefaultTesterConfiguration("");
  
  // third argument: path of the parameter set for geometry configuration
  // (optional; default: "services.Geometry" from the inherited object)
  if (++iParam < argc) config.SetGeometryParameterSetPath(argv[iParam]);
  
  //
  // testing environment setup
  //
  using namespace lar::standalone;
  SBNDGeometryTestEnvironment TestEnvironment(config);
  auto const wireReadout = SetupReadout<geo::WireReadoutSorterSBND>(TestEnvironment.ServiceParameters("WireReadout"),
                                                                    TestEnvironment.Geometry());
  auto const auxDetGeom = SetupAuxDetGeometry(TestEnvironment.ServiceParameters("AuxDetGeometry"));

  //
  // run the test algorithm
  //
  
  // 1. we initialize it from the environment,
  geo::GeometryTestAlg Tester(TestEnvironment.Geometry(),
                              wireReadout.get(),
                              auxDetGeom.get(),
                              TestEnvironment.TesterParameters());
  
  // 2. then we run it!
  unsigned int nErrors = Tester.Run();
  
  // 3. And finally we cross fingers.
  if (nErrors > 0) {
    mf::LogError("geometry_test_SBND") << nErrors << " errors detected!";
  }
  
  return nErrors;
} // main()
