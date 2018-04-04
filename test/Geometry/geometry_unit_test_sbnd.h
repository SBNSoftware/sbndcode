/**
 * @file   geometry_unit_test_sbnd.h
 * @brief  Class for objects initializing SBND geometry
 * @date   June 6, 2016
 * @author petrillo@fnal.gov
 * 
 * Provides an environment for easy set up of SBND-aware tests.
 * Keep in mind that the channel mapping algorithm must be hard-coded and, if
 * using Boost unit test, the configuration file location must be hard coded too
 * (or you can use the provided configuration).
 * 
 * For an example of usage, see larcore/test/Geometry/geometry_test.cxx
 */

#ifndef TEST_GEOMETRY_UNIT_TEST_SBND_H
#define TEST_GEOMETRY_UNIT_TEST_SBND_H

// LArSoft libraries
#include "test/Geometry/geometry_unit_test_base.h"
// #include "sbndcode/Geometry/ChannelMapSBNDAlg.h"
#include "larcorealg/Geometry/ChannelMapStandardAlg.h"

// C/C++ standard libraries
#include <string>


namespace geo {
  class ChannelMapStandardAlg;
} // namespace geo

namespace sbnd {
  
  /// Namespace including SBND-specific testing
  namespace testing {
    
    /** ************************************************************************
     * @brief Class holding the configuration for a SBND fixture
     * @tparam CHANNELMAP the class used for channel mapping
     * @see BasicGeometryEnvironmentConfiguration
     *
     * This class needs to be fully constructed by the default constructor
     * in order to be useful as Boost unit test fixture.
     * It is supposed to be passed as a template parameter to another class
     * that can store an instance of it and extract configuration information
     * from it.
     * 
     * This class should be used with ChannelMapStandardAlg.
     * 
     * We reuse BasicGeometryEnvironmentConfiguration as base class and then we
     * fix its setup.
     */
    template <typename CHANNELMAP = geo::ChannelMapStandardAlg>
    struct SBNDGeometryEnvironmentConfiguration:
      public ::testing::BasicGeometryEnvironmentConfiguration<CHANNELMAP>
    {
      // remember that BasicGeometryEnvironmentConfiguration is not polymorphic
      using base_t
        = ::testing::BasicGeometryEnvironmentConfiguration<CHANNELMAP>;
      
      /// Default constructor
      SBNDGeometryEnvironmentConfiguration() { SBNDdefaultInit(); }
      
      /// Constructor; accepts the name as parameter
      SBNDGeometryEnvironmentConfiguration(std::string name):
        SBNDGeometryEnvironmentConfiguration()
        { base_t::SetApplicationName(name); }
      
      
        private:
      void SBNDdefaultInit()
        {
          // overwrite the configuration that happened in the base class:
          base_t::SetApplicationName("SBNDGeometryTest");
          base_t::SetDefaultGeometryConfiguration(R"(
            SurfaceY: 130e2 # in cm, vertical distance to the surface
            Name:     "sbndv1"
            GDML:     "sbnd_v01_00.gdml"
            ROOT:     "sbnd_v01_00.gdml"
            SortingParameters: {}
            )");
        }
    }; // class SBNDGeometryEnvironmentConfiguration<>
    
    
  } // namespace testing
} // namespace sbnd

#endif // TEST_GEOMETRY_UNIT_TEST_SBND_H
