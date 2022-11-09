
set (
          MODULE_LIBRARIES larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation lardataobj_Simulation
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           lardata_Utilities
                           larevt_Filters
                           lardataobj_RawData
                           lardataobj_RecoBase
                           lardataobj_AnalysisBase
                           lardata_RecoObjects
                           larpandora_LArPandoraInterface
                           nusimdata::SimulationBase
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art_root_io::tfile_support ROOT::Core
                           art_root_io::TFileService_service
                           art::Persistency_Common canvas
                           art::Persistency_Provenance canvas
                           art::Utilities canvas
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
                           sbnobj::Common_CRT
                           sbndcode_CRTUtils
                           sbndcode_GeoWrappers
        ) 

cet_build_plugin(CRTT0Matching art::module SOURCE CRTT0Matching_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(CRTTrackMatching art::module SOURCE CRTTrackMatching_module.cc LIBRARIES ${MODULE_LIBRARIES})


install_headers()
install_fhicl()
install_source()

