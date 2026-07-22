
cet_build_plugin( SBNDCRTGeometryTest art::module LIBRARIES
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Framework_Services_Optional_RandomNumberGenerator_service
  art::Persistency_Provenance
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
  NO_INSTALL
)

cet_transitive_paths(FHICL_DIR BINARY IN_TREE)
cet_test_env_prepend(FHICL_FILE_PATH "." ${TRANSITIVE_PATHS_WITH_FHICL_DIR} )
cet_transitive_paths(LIBRARY_DIR BINARY IN_TREE)
cet_test_env_prepend(CET_PLUGIN_PATH ${TRANSITIVE_PATHS_WITH_LIBRARY_DIR})
cet_transitive_paths(GDML_DIR BINARY IN_TREE)
cet_test_env_prepend(FW_SEARCH_PATH ${TRANSITIVE_PATHS_WITH_GDML_DIR})

cet_test(SBNDCRTGeometryTest_1 HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config SBNDCRTGeometryTest.fcl -n 1
  DATAFILES
    SBNDCRTGeometryTest.fcl
)
