/**
 * @file   geometry_unit_test_lar1nd.h
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

#ifndef TEST_GEOMETRY_UNIT_TEST_LAR1ND_H
#define TEST_GEOMETRY_UNIT_TEST_LAR1ND_H

// LArSoft libraries
#include "test/Geometry/geometry_unit_test_base.h"
#include "lar1ndcode/Geo/ChannelMaplar1ndAlg.h"

// C/C++ standard libraries
#include <string>


namespace geo {
  class ChannelMapStandardAlg;
} // namespace geo

namespace lar1nd {
  
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
    template <typename CHANNELMAP = geo::ChannelMaplar1ndAlg>
    struct LAr1NDGeometryEnvironmentConfiguration:
      public ::testing::BasicGeometryEnvironmentConfiguration<CHANNELMAP>
    {
      // remember that BasicGeometryEnvironmentConfiguration is not polymorphic
      using base_t
        = ::testing::BasicGeometryEnvironmentConfiguration<CHANNELMAP>;
      
      /// Default constructor
      LAr1NDGeometryEnvironmentConfiguration() { LAr1NDdefaultInit(); }
      
      /// Constructor; accepts the name as parameter
      LAr1NDGeometryEnvironmentConfiguration(std::string name):
        LAr1NDGeometryEnvironmentConfiguration()
        { base_t::SetApplicationName(name); }
      
      
        private:
      void LAr1NDdefaultInit()
        {
          // overwrite the configuration that happened in the base class:
          base_t::SetApplicationName("LAr1NDGeometryTest");
          base_t::SetDefaultGeometryConfiguration(R"(
            SurfaceY: 130e2 # in cm, vertical distance to the surface
            Name:     "lar1ndv0"
            GDML:     "lar1ndv0.gdml"
            ROOT:     "lar1ndv0.gdml"
            SortingParameters: {}
            )");
        }
    }; // class LAr1NDGeometryEnvironmentConfiguration<>
    
    
  } // namespace testing
} // namespace lar1nd

#endif // TEST_GEOMETRY_UNIT_TEST_LAR1ND_H
