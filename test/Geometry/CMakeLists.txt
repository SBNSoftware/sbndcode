
# geometry test on SBND geometries

# these tests run a FCL file and fail only if lar exits with a bad exit code;
# they are fast and always run.
# FCL files need to be copied to the test area (DATAFILES directive) since they
# are not installed.
cet_test(geometry_lar1nd HANDBUILT
  DATAFILES test_geometry_lar1nd.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./test_geometry_lar1nd.fcl
)


#
# the following are tests equivalent to the module ones, but with no framework
#
# unit test (use the hard-coded configuration for SBND v0 geometry)
cet_test(geometry_lar1nd_test
  SOURCES geometry_lar1nd_test.cxx
  DATAFILES test_geometry_lar1nd.fcl
  TEST_ARGS ./test_geometry_lar1nd.fcl
  LIBRARIES larcore_Geometry
            GeometryTestLib
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
	    ${ROOT_CORE}
)


# unit test (use the hard-coded configuration for SBND v0 geometry)
# this uses BOOST for the test
cet_test(geometry_iterator_lar1nd_test
  SOURCES geometry_iterator_lar1nd_test.cxx
  LIBRARIES larcore_Geometry
            GeometryTestLib
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
	    ${ROOT_CORE}
  USE_BOOST_UNIT
)


# unit test (use the hard-coded configuration for SBND v0 geometry)
cet_test(geometry_iterator_loop_lar1nd_test
  SOURCES geometry_iterator_loop_lar1nd_test.cxx
  LIBRARIES larcore_Geometry
            GeometryTestLib
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
	    ${ROOT_CORE}
)


# FCL files need to be copied to the test area (DATAFILES directive) since they
# are not installed.
cet_test(dump_channel_map_lar1nd_test HANDBUILT
  DATAFILES dump_lar1nd_channelmap.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./dump_lar1nd_channelmap.fcl
)

