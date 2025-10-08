cet_enable_asserts()

set(    MODULE_LIBRARIES
  	art_root_io::TFileService_service
      art_root_io::tfile_support
)
set(	LIB_LIBRARIES
      art::Framework_Services_Registry
      messagefacility::MF_MessageLogger
      lardata::Utilities
	larevt::CalibrationDBI_IOVData
	larevt::CalibrationDBI_Providers
)
set(	SERVICE_LIBRARIES
      sbndcode_Calibration_PDSDatabaseInterface
      larcore::Geometry_Geometry_service
      lardata::DetectorClocksService
)

file(GLOB lib_srcs *.cxx)

art_make_library( SOURCE ${lib_srcs} LIBRARIES PUBLIC ${LIB_LIBRARIES})

cet_build_plugin( PMTCalibrationDatabaseService art::service LIBRARIES PUBLIC ${SERVICE_LIBRARIES})

install_headers()
install_fhicl()
install_source()