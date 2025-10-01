cet_enable_asserts()

set(    LIB_LIBRARIES
      art::Framework_Services_Registry
      messagefacility::MF_MessageLogger
      lardata::Utilities
        larevt::CalibrationDBI_IOVData
        larevt::CalibrationDBI_Providers
)
set(    SERVICE_LIBRARIES
      sbnobj::SBND_CRT
      sbndcode_Calibration_CRTDatabaseInterface
      larcore::Geometry_Geometry_service
      lardata::DetectorClocksService
)

file(GLOB lib_srcs *.cxx)

art_make_library( SOURCE ${lib_srcs} LIBRARIES PUBLIC ${LIB_LIBRARIES})

cet_build_plugin( CRTCalibrationDatabaseService art::service LIBRARIES PUBLIC ${SERVICE_LIBRARIES})

install_headers()
install_fhicl()
install_source()
