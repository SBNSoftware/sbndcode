include(CetTest)
# geometry test on SBND geometries
cet_transitive_paths(FHICL_DIR BINARY IN_TREE)
cet_test_env_prepend(FHICL_FILE_PATH "." ${TRANSITIVE_PATHS_WITH_FHICL_DIR} ${CMAKE_CURRENT_LIST_DIR})
cet_transitive_paths(LIBRARY_DIR BINARY IN_TREE)
cet_test_env_prepend(CET_PLUGIN_PATH ${TRANSITIVE_PATHS_WITH_LIBRARY_DIR})
cet_transitive_paths(GDML_DIR BINARY IN_TREE)
cet_test_env_prepend(FW_SEARCH_PATH ${TRANSITIVE_PATHS_WITH_GDML_DIR})
# these tests run a FCL file and fail only if lar exits with a bad exit code;
# they are fast and always run.
# FCL files need to be copied to the test area (DATAFILES directive) since they
# are not installed.
cet_test(geometry_sbnd HANDBUILT
  DATAFILES test_geometry_sbnd.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./test_geometry_sbnd.fcl
#  OPTIONAL_GROUPS Broken
)


#
# the following are tests equivalent to the module ones, but with no framework
#
# unit test (uses configuration in test_geometry_sbnd.fcl)
cet_test(geometry_sbnd_test
  SOURCE geometry_sbnd_test.cxx
  DATAFILES test_geometry_sbnd.fcl
  TEST_ARGS ./test_geometry_sbnd.fcl
  LIBRARIES sbndcode_Geometry
            larcorealg::Geometry
            larcorealg::GeometryTestLib
            messagefacility::MF_MessageLogger
            fhiclcpp::fhiclcpp
            cetlib_except::cetlib_except
	    ROOT::Core
#  OPTIONAL_GROUPS Broken
)


# unit test (uses configuration in test_geometry_sbnd.fcl)
# this uses BOOST for the test
cet_test(geometry_iterator_sbnd_test
  DATAFILES test_geometry_sbnd.fcl
  TEST_ARGS ./test_geometry_sbnd.fcl
  SOURCE geometry_iterator_sbnd_test.cxx
  LIBRARIES sbndcode_Geometry
            larcorealg::Geometry
            larcorealg::GeometryTestLib
            messagefacility::MF_MessageLogger
            fhiclcpp::fhiclcpp
            cetlib_except::cetlib_except
	    ROOT::Core
  USE_BOOST_UNIT
)


# unit test (uses configuration in test_geometry_sbnd.fcl)
cet_test(geometry_iterator_loop_sbnd_test
  SOURCE geometry_iterator_loop_sbnd_test.cxx
  DATAFILES test_geometry_sbnd.fcl
  TEST_ARGS ./test_geometry_sbnd.fcl
  LIBRARIES sbndcode_Geometry
            larcorealg::Geometry
            larcorealg::GeometryTestLib
            messagefacility::MF_MessageLogger
            fhiclcpp::fhiclcpp
            cetlib_except::cetlib_except
	    ROOT::Core
)


# FCL files need to be copied to the test area (DATAFILES directive) since they
# are not installed.
cet_test(dump_channel_map_sbnd_test HANDBUILT
  DATAFILES dump_sbnd_channelmap.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./dump_sbnd_channelmap.fcl
)


install_headers()
add_subdirectory(CRTGeometry)
