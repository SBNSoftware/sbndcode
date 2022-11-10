add_subdirectory(alg)

set(     MODULE_LIBRARIES
            sbndcode_MCTruthExtractor_alg
            larcore::Geometry_Geometry_service
            art_root_io::TFileService_service
            art::Framework_Services_Optional_RandomNumberGenerator_service
            larcorealg::Geometry
            larsim::Simulation
            nug4::ParticleNavigation
            lardataobj::Simulation
            lardata::Utilities
            larevt::Filters
            lardataobj::RawData
            nugen::NuReweight_art
            nusimdata::SimulationBase
            art::Framework_Core
            art::Framework_Principal
            art::Framework_Services_Registry
            art_root_io::tfile_support
            ROOT::Core
            ROOT::Tree
            art::Persistency_Common
            art::Persistency_Provenance
            art::Utilities
            canvas::canvas
            messagefacility::MF_MessageLogger
            fhiclcpp::fhiclcpp
            cetlib::cetlib
            CLHEP::CLHEP
          ) 
cet_build_plugin( NuAna art::module SOURCE NuAna_module.cc LIBRARIES ${MODULE_LIBRARIES})
target_include_directories (sbndcode_MCTruthExtractor_NuAna_module PRIVATE $ENV{GENIE_INC}/GENIE )
install_headers()
install_fhicl()
