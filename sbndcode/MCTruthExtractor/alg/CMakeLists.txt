include_directories ( $ENV{GENIE_INC}/GENIE )

art_make(
            LIB_LIBRARIES
            larcore_Geometry_Geometry_service
            larcorealg_Geometry
            nugen::NuReweight_art
            nugen::NuReweight
            nusimdata::SimulationBase
            art::Framework_Core
            art::Framework_Principal
            art::Framework_Services_Registry
            art_root_io::tfile_support ROOT::Core
            art::Persistency_Common
            art::Persistency_Provenance
            art::Utilities
            canvas::canvas
            messagefacility::MF_MessageLogger

            fhiclcpp::fhiclcpp
            cetlib::cetlib cetlib_except
            CLHEP::CLHEP
            
            ${ROOT_BASIC_LIB_LIST}
          )

install_headers()

