art_make_library( SOURCE NuAnaAlg.cxx
            LIBRARIES
            larcore::Geometry_Geometry_service
            larcorealg::Geometry
            nugen::NuReweight_art
            nugen::NuReweight
            nugen::EventGeneratorBase_GENIE
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
            ROOT::Core
          ) 

install_headers()

