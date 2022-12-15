
set (
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

cet_build_plugin(CRTT0Matching art::module SOURCE CRTT0Matching_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(CRTTrackMatching art::module SOURCE CRTTrackMatching_module.cc LIBRARIES ${MODULE_LIBRARIES})


install_headers()
install_fhicl()
install_source()

