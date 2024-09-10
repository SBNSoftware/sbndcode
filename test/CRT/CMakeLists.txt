
cet_test(crt_detsimalg_test
  SOURCE crt_detsimalg_test.cxx
  DATAFILES test_crt_detsimalg.fcl
  TEST_ARGS ./test_crt_detsimalg.fcl
  LIBRARIES sbndcode_CRT_CRTSimulation
            sbndcode_Geometry
            larcorealg::Geometry
            larcorealg::GeometryTestLib
            messagefacility::MF_MessageLogger
            fhiclcpp::fhiclcpp
            cetlib_except::cetlib_except
      ROOT::Core
#  OPTIONAL_GROUPS Broken
)

install_headers()
