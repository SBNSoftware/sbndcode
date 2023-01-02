add_subdirectory(Dumpers)

cet_enable_asserts()

art_make_library(SUBDIRS details

        LIBRARIES
                lardataalg::DetectorInfo
                sbnobj::Common_PMT_Data
                messagefacility::MF_MessageLogger
                fhiclcpp::fhiclcpp
                cetlib_except::cetlib_except
                cetlib::cetlib
                ROOT::RIO
                ROOT::Core
                art::Framework_Services_Registry
                FFTW3::FFTW3
)

set(        TOOL_LIBRARIES
                sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
                artdaq_core::artdaq-core_Utilities
                sbndcode::Decoders_DecoderTools
                sbndcode::Decoders_DecoderTools_Dumpers
                sbndcode::Decoders_DecoderTools_details
                sbndcode::Utilities
                larcorealg::Geometry
                larcore::Geometry_Geometry_service
                lardata::Utilities
                lardataalg::DetectorInfo
                larevt::Filters
                lardataobj::RawData
                lardata::ArtDataHelper
                ROOT::Core
                wda::wda
                art::Framework_Core
                art::Framework_Principal
                art::Framework_Services_Registry
                art_root_io::tfile_support
                art_root_io::TFileService_service
                art::Utilities
                art::Persistency_Common canvas::canvas
                art::Persistency_Provenance
                art::Utilities
                messagefacility::MF_MessageLogger
                messagefacility::headers
                fhiclcpp::fhiclcpp
                cetlib::cetlib cetlib_except::cetlib_except
                Boost::program_options
        )

install_headers()
install_fhicl()
install_source()
