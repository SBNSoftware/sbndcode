cet_enable_asserts()

set(    MODULE_LIBRARIES
        lardataobj::RecoBase
	    sbnobj::Common_Trigger
        larcore::Geometry_Geometry_service
        lardata::DetectorClocksService
	    art_root_io::TFileService_service
        art_root_io::tfile_support
        lardataobj::RawData
        ROOT::Tree
)
set(	LIB_LIBRARIES
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        lardataobj::RecoBase
        lardata::Utilities
        sbnobj::Common_PMT_Data
	    larevt::CalibrationDBI_IOVData
	    larevt::CalibrationDBI_Providers
)
set(	SERVICE_LIBRARIES
	    sbndcode_DatabaseInterface
        larcore::Geometry_Geometry_service
        lardata::DetectorClocksService
)

file(GLOB lib_srcs *.cxx)

art_make_library( SOURCE ${lib_srcs} LIBRARIES PUBLIC ${LIB_LIBRARIES})

cet_build_plugin( PMTCalibrationDatabaseService art::service LIBRARIES PUBLIC ${SERVICE_LIBRARIES})

install_headers()
install_fhicl()
install_source()