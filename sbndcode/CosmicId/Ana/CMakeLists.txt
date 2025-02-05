
set(
          MODULE_LIBRARIES larcorealg::Geometry
                           larcore::Geometry_Geometry_service
                           larsim::Simulation
                           lardataobj::Simulation
                           larsim::MCCheater_BackTrackerService_service
                           larsim::MCCheater_ParticleInventoryService_service
                           lardata::Utilities
                           larevt::Filters
                           lardataobj::RawData
                           lardataobj::RecoBase
                           lardataobj::AnalysisBase
                           lardata::RecoObjects
                           larpandora::LArPandoraInterface
                           nusimdata::SimulationBase
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art_root_io::tfile_support
                           art_root_io::TFileService_service
                           art::Persistency_Common
                           art::Persistency_Provenance
                           art::Utilities
                           canvas::canvas
                           messagefacility::MF_MessageLogger
                           fhiclcpp::fhiclcpp
                           ROOT::Geom
                           ROOT::XMLIO
                           ROOT::Gdml
                           ROOT::Core
                           ROOT::Spectrum
                           ROOT::RooFit
                           ROOT::RooFitCore
                           sbndcode_RecoUtils
                           sbndcode_CRT
                           sbndcode_CRTUtils
                           sbndcode_GeoWrappers
                           sbndcode_CosmicIdAlgs
                           sbndcode_CosmicIdUtils
                           larreco::RecoAlg
        )

cet_build_plugin(CosmicIdAna art::module SOURCE CosmicIdAna_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(CosmicIdTree art::module SOURCE CosmicIdTree_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(CRTHitCosmicIdAna art::module SOURCE CRTHitCosmicIdAna_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(CRTTrackCosmicIdAna art::module SOURCE CRTTrackCosmicIdAna_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(StoppingCosmicIdAna art::module SOURCE StoppingCosmicIdAna_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()

