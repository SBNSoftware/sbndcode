cet_enable_asserts()

set(    MODULE_LIBRARIES
        lardataobj::RecoBase
        lardataobj::RawData
        larcore::Geometry_Geometry_service
        lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
	lardataalg::DetectorInfo
)
set(    LIB_LIBRARIES
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        lardataobj::RecoBase
        lardataobj::RawData
        lardata::Utilities
        wda::wda
        sbnobj::Common_PMT_Data
	lardataalg::DetectorInfo
)

set(    SERVICE_LIBRARIES
        larcore::Geometry_Geometry_service
        lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
)

file(GLOB lib_srcs *.cxx)

art_make_library(
	SOURCE
		"PMTDecoderUtils.cxx"
                ${lib_srcs} LIBRARIES PUBLIC ${LIB_LIBRARIES}
	LIBRARIES
		artdaq_core::artdaq-core_Data
		cetlib_except::cetlib_except
                art_root_io::RootDB
                SQLite::SQLite3
                lardata::Utilities
                canvas::canvas
                cetlib::cetlib
                fhiclcpp::fhiclcpp
                hep_concurrency::hep_concurrency # required by fhicl::ParameterSetRegistry
                ROOT::Tree
                ROOT::Core
                ROOT::RIO
	)

install_headers()
install_source()
install_fhicl()
