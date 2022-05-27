
simple_plugin( SBNDCRTGeometryTest "module"
  larcorealg_Geometry
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Framework_Services_Optional_RandomNumberGenerator_service
  art::Persistency_Provenance
  canvas::canvas
  messagefacility::MF_MessageLogger

  fhiclcpp::fhiclcpp
  cetlib::cetlib cetlib_except
  NO_INSTALL
)

cet_test(SBNDCRTGeometryTest_1 HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config SBNDCRTGeometryTest.fcl -n 1
  DATAFILES
    SBNDCRTGeometryTest.fcl
)
