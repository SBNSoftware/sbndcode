/**
 * @file   geometry_iterator_lar1nd_test.cxx
 * @brief  Unit test for geometry iterators on SBND detector
 * @date   June 6th, 2016
 * @author petrillo@fnal.gov
 * 
 * Usage: just run the executable.
 * Boost unit testing environment keeps the arguments secret anyway.
 */

// Boost test libraries; defining this symbol tells boost somehow to generate
// a main() function; Boost is pulled in by geometry_boost_unit_test_base.h
#define BOOST_TEST_MODULE GeometryIteratorTestLAr1ND

// LArSoft libraries
#include "test/Geometry/geometry_unit_test_lar1nd.h"
#include "test/Geometry/GeometryIteratorTestAlg.h"
#include "test/Geometry/boost_unit_test_base.h"
#include "larcore/Geometry/ChannelMapStandardAlg.h"

//------------------------------------------------------------------------------
//---  The test environment
//---

// we define here all the configuration that is needed;
// in the specific, the type of the channel mapping and a proper test name,
// used for output only; we use LAr1NDGeometryEnvironmentConfiguration
// as base class, that is already configured to use SBND geometry.
// We wrap it in testing::BoostCommandLineConfiguration<> so that it can learn
// the configuration file name from the command line.
struct LAr1NDGeometryConfiguration:
  public testing::BoostCommandLineConfiguration<
    lar1nd::testing::LAr1NDGeometryEnvironmentConfiguration
      <geo::ChannelMapStandardAlg>
    >
{
  /// Constructor: overrides the application name; ignores command line
  LAr1NDGeometryConfiguration()
    { SetApplicationName("GeometryIteratorUnitTest"); }
}; // class LAr1NDGeometryConfiguration

/*
 * Our fixture is based on GeometryTesterFixture, configured with the object
 * above.
 * It provides:
 * - `Tester`, a configured instance of the test algorithm.
 */
class LAr1NDGeometryIteratorTestFixture:
  private testing::GeometryTesterEnvironment<LAr1NDGeometryConfiguration>
{
    public:
  geo::GeometryIteratorTestAlg Tester;
  
  /// Constructor: initialize the tester with the Geometry from base class
  LAr1NDGeometryIteratorTestFixture(): Tester(TesterParameters())
    { Tester.Setup(*Geometry()); }

}; // class LAr1NDGeometryIteratorTestFixture



//------------------------------------------------------------------------------
//---  The tests
//---

BOOST_FIXTURE_TEST_SUITE
  (GeometryIteratorsLAr1ND, LAr1NDGeometryIteratorTestFixture)
// BOOST_GLOBAL_FIXTURE(LAr1NDGeometryIteratorTestFixture)


BOOST_AUTO_TEST_CASE( AllTests )
{
  Tester.Run();
} // BOOST_AUTO_TEST_CASE( AllTests )

/*
BOOST_AUTO_TEST_CASE( CryostatIteratorsTest )
{
  Tester.CryostatIteratorsTest();
} // BOOST_AUTO_TEST_CASE( CryostatIteratorsTest )



BOOST_AUTO_TEST_CASE( TPCIteratorsTest )
{
  Tester.TPCIteratorsTest();
} // BOOST_AUTO_TEST_CASE( TPCIteratorsTest )



BOOST_AUTO_TEST_CASE( PlaneIteratorsTest )
{
  Tester.PlaneIteratorsTest();
} // BOOST_AUTO_TEST_CASE( PlaneIteratorsTest )



BOOST_AUTO_TEST_CASE( WireIteratorsTest )
{
  Tester.WireIteratorsTest();
} // BOOST_AUTO_TEST_CASE( WireIteratorsTest )
*/

BOOST_AUTO_TEST_SUITE_END()

