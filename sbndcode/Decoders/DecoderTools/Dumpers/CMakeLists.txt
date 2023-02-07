cet_enable_asserts()

art_dictionary(DICTIONARY_LIBRARIES
                 lardataobj::RecoBase
                 lardataalg::headers
                 canvas::canvas
              )

set(    MODULE_LIBRARIES
        lardataobj::RecoBase
        lardataobj::RawData
        larcore::Geometry_Geometry_service
        lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
	lardataalg::DetectorInfo
	wda::wda
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
        wda::wda
)

file(GLOB lib_srcs *.cxx)

art_make_library(
	SOURCE
		"FragmentDumper.cxx"
                ${lib_srcs} LIBRARIES PUBLIC ${LIB_LIBRARIES}
	LIBRARIES
                art_root_io::RootDB
                art::Framework_Services_Registry
		artdaq_core::artdaq-core_Data
                canvas::canvas
                cetlib::cetlib
		cetlib_except::cetlib_except
                fhiclcpp::fhiclcpp
                FFTW3::FFTW3
                hep_concurrency::hep_concurrency # required by fhicl::ParameterSetRegistry
                lardata::Utilities
                lardataalg::DetectorInfo
                messagefacility::MF_MessageLogger
                ROOT::Tree
                ROOT::Core
                ROOT::RIO
                sbnobj::Common_PMT_Data
                SQLite::SQLite3

	)

install_headers()
install_source()
install_fhicl()
