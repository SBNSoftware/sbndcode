
set(      MODULE_LIBRARIES larcorealg::Geometry
                           larcore::Geometry_Geometry_service
                           larsim::Simulation lardataobj::Simulation
                           larsim::MCCheater_BackTrackerService_service
                           larsim::MCCheater_ParticleInventoryService_service
                           lardata::Utilities
                           larevt::Filters
                           lardataobj::RawData
                           lardataobj::RecoBase
			   larreco::Calorimetry
                           larreco::RecoAlg
                           lardata::RecoObjects
                           larpandora::LArPandoraInterface
                           sbndcode_CRTUtils
                           sbnobj::Common_CRT
                           nusimdata::SimulationBase
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art_root_io::tfile_support
                           art_root_io::TFileService_service
                           art::Persistency_Common
                           art::Persistency_Provenance
                           art::Utilities
                           messagefacility::MF_MessageLogger
                           ROOT::Core
                           ROOT::Tree
                           fhiclcpp::fhiclcpp
                           ROOT::Geom
                           ROOT::XMLIO
                           ROOT::Gdml
                           sbndcode_RecoUtils
                           sbndcode_OpDetSim
        )

cet_build_plugin ( TPCRecohitAna art::module
  larcore::Geometry_Geometry_service
  lardataobj::Simulation
  lardataobj::RawData
  lardataobj::RecoBase
  ROOT::Core
  ROOT::Tree
  art::Framework_Core
  art::Framework_Principal
  art::Utilities
  art_root_io::tfile_support
  art_root_io::TFileService_service
  canvas::canvas
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
  messagefacility::MF_MessageLogger

  BASENAME_ONLY
              )


cet_build_plugin(AnalysisTree art::Module SOURCE AnalysisTree_module.cc LIBRARIES ${MODULE_LIBRARIES} )
# install_headers()
install_fhicl()
install_source()

