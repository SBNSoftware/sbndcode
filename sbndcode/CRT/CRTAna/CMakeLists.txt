
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
                           art_root_io::tfile_support ROOT::Core
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
                           ROOT::Tree
                           ROOT::Spectrum
                           ROOT::RooFit
                           ROOT::RooFitCore
                           sbndcode_RecoUtils
                           sbndcode_CRT
                           sbnobj::Common_CRT
                           sbndcode_CRTUtils
                           sbndcode_GeoWrappers
        )

cet_build_plugin(CRTDetSimAna art::module SOURCE CRTDetSimAna_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(CRTFullRecoAna art::module SOURCE CRTFullRecoAna_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(CRTHitRecoAna art::module SOURCE CRTHitRecoAna_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(CRTT0MatchingAna art::module SOURCE CRTT0MatchingAna_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(CRTTrackMatchingAna art::module SOURCE CRTTrackMatchingAna_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(CRTTrackRecoAna art::module SOURCE CRTTrackRecoAna_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()

