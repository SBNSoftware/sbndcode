art_make_library(    LIBRARY_NAME sbndcode_CRTUtils
                   SOURCE
                           CRTBackTracker.cc
                           CRTCommonUtils.cc
                           CRTEventDisplay.cc
                           CRTHitRecoAlg.cc
                           CRTT0MatchAlg.cc
                           CRTTrackMatchAlg.cc
                           CRTTrackRecoAlg.cc
                           TPCGeoUtil.cc
                 LIBRARIES larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation lardataobj_Simulation
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           lardata_Utilities
                           larevt_Filters
                           lardataobj_RawData
                           lardataobj_RecoBase
                           lardata_RecoObjects  
                           lardataalg_DetectorInfo
                           larpandora_LArPandoraInterface
                           nusimdata::SimulationBase
                           cetlib::cetlib cetlib_except
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
                           sbndcode_CRT
                           sbnobj::Common_CRT
                           sbnobj::SBND_CRT
                           sbndcode_GeoWrappers
                           sbndcode_RecoUtils
        ) 

install_headers()
install_fhicl()
install_source()
