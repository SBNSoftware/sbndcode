cet_enable_asserts()

set(          MODULE_LIBRARIES
			sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
                        artdaq_core::artdaq-core_Utilities
                        larcorealg::Geometry
                        lardataobj::RawData
                        larcore::Geometry_Geometry_service
                        lardata::Utilities
                        larevt::Filters
                        lardataobj::RawData
                        lardata::ArtDataHelper
)

cet_build_plugin(DaqDecoderSBNDPMT art::module LIBRARIES
	sbndcode::Decoders_DecoderTools_Dumpers
	sbndcode::Decoders_DecoderTools
	sbnobj::Common_PMT_Data
	sbnobj::Common_Trigger
	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
	lardataobj::RawData
	larcore::Geometry_Geometry_service
	artdaq_core::artdaq-core_Data
	art_root_io::TFileService_service
	art_root_io::tfile_support
	art::Framework_Services_Registry
	messagefacility::MF_MessageLogger
	fhiclcpp::fhiclcpp
	cetlib_except::cetlib_except
	ROOT::Tree
	lardata::Utilities
	)

cet_build_plugin(PMTconfigurationExtraction art::module LIBRARIES
	sbndcode::Decoders_DecoderTools_Dumpers
	sbnobj::Common_PMT_Data
	sbnobj::Common_Trigger
	art::Framework_Services_Registry
	fhiclcpp::fhiclcpp
	cetlib_except::cetlib_except
	)

cet_build_plugin(DumpPMTconfiguration art::module LIBRARIES
	sbnobj::Common_PMT_Data
        lardataobj::RawData
	messagefacility::MF_MessageLogger
	fhiclcpp::fhiclcpp
	cetlib_except::cetlib_except
	)

install_headers()
install_fhicl()
install_source()

# Add our tools directory
add_subdirectory(DecoderTools)

