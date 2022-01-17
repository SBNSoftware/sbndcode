add_subdirectory(alg)
include_directories ( $ENV{GENIE_INC}/GENIE )

art_make(
            LIBRARY_NAME Ntuples
            MODULE_LIBRARIES
            sbndcode_MCTruthExtractor_alg

            larcore_Geometry_Geometry_service
            # TimeService_service
            art_root_io::TFileService_service
            art::Framework_Services_Optional_RandomNumberGenerator_service
            larcorealg_Geometry
            larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
            lardata_Utilities
            larevt_Filters
            lardataobj_RawData
            nugen::NuReweight_art
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
install_fhicl()
