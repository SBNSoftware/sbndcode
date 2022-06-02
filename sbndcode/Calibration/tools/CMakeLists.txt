cet_enable_asserts()

include_directories($ENV{EIGEN_INC})

art_make( TOOL_LIBRARIES lardataobj_RecoBase
                         lardataobj_Simulation
                         larcorealg_Geometry
                         larcore_Geometry_Geometry_service
                         lardata_RecoObjects
                         lardata_Utilities
                         nusimdata::SimulationBase
                         larsim_Simulation
                         nusimdata::SimulationBase
                         nuevdb::EventDisplayBase
                         art::Framework_Core
                         art::Framework_Principal
                         art::Framework_Services_Registry
                         ${ART_FRAMEWORK_SERVICES_BASIC}
                         art_root_io::tfile_support ROOT::Core
                         art_root_io::TFileService_service
                         art::Framework_Services_Optional_RandomNumberGenerator_service
                         art::Persistency_Common
                         art::Persistency_Provenance
                         art::Utilities
                         canvas::canvas
                         messagefacility::MF_MessageLogger
                         
                         fhiclcpp::fhiclcpp
                         cetlib::cetlib
                         cetlib_except::cetlib_except
                         Boost::filesystem
                         
                         CLHEP::CLHEP
                         ${ROOT_BASIC_LIB_LIST}
                         ${ROOT_MINUIT}
        )

install_headers()
install_source()
install_fhicl()
