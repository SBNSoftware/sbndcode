cet_enable_asserts()

set     ( TOOL_LIBRARIES lardataobj::RecoBase
                         lardataobj::Simulation
                         larcorealg::Geometry
                         larcore::Geometry_Geometry_service
                         lardata::RecoObjects
                         lardata::Utilities
                         larevt::CalibrationDBI_Providers
                         nusimdata::SimulationBase
                         larsim::Simulation
                         nusimdata::SimulationBase
                         nuevdb::EventDisplayBase
                         art::Framework_Core
                         art::Framework_Principal
                         art::Framework_Services_Registry
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
                         ROOT::Core
                         ROOT::Minuit
                         Eigen3::Eigen
        )
cet_build_plugin(TrackHitEfficiencyAnalysis art::tool SOURCE TrackHitEfficiencyAnalysis_tool.cc LIBRARIES ${TOOL_LIBRARIES})

install_headers()
install_source()
install_fhicl()
