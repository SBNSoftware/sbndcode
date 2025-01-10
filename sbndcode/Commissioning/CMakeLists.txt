add_subdirectory(fcls)

set(
   MODULE_LIBRARIES
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
        larsim::Simulation
        lardataobj::Simulation
        larsim::MCCheater_BackTrackerService_service
        larsim::MCCheater_ParticleInventoryService_service
        lardata::Utilities
        larevt::Filters
        lardataobj::RawData
        lardataobj::RecoBase
        lardataobj::MCBase
        lardataobj::Simulation
        larreco::RecoAlg
        lardata::RecoObjects
        larpandora::LArPandoraInterface
        sbndcode_CRTUtils
        sbnobj::Common_CRT
        sbnobj::SBND_CRT
        sbndaq_artdaq_core::sbndaq-artdaq-core_Obj_SBND
        nusimdata::SimulationBase
        art::Framework_Core
        art::Framework_Principal
        art::Framework_Services_Registry
        art_root_io::tfile_support ROOT::Core
        art_root_io::TFileService_service
        art::Persistency_Common
        art::Persistency_Provenance
        art::Utilities
        messagefacility::MF_MessageLogger
        ROOT::Core
        fhiclcpp::fhiclcpp
        ROOT::Geom
        ROOT::XMLIO
        ROOT::Gdml
        sbndcode_RecoUtils
        sbndcode_OpDetSim
        sbndcode_GeoWrappers
        sbndcode_ChannelMaps_TPC_TPCChannelMapService_service
        sbndcode_Utilities_DigitalNoiseChannelStatus_service
        hep_hpc_hdf5
        ${HDF5_LIBRARIES}
)



set (
  Ringing_LIBRARIES
                    sbndcode_OpDetSim
                    larcore::Geometry_Geometry_service
                    lardataobj::Simulation
                    lardata::Utilities
                    lardataobj::RawData
                    lardata::Utilities_LArFFT_service
                    lardataobj::RecoBase
                    lardata::Utilities_LArFFT_service
                    lardataobj::AnalysisBase
                    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
                    sbndcode_Utilities_SignalShapingServiceSBND_service
                    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
                    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
                    sbndaq_artdaq_core::sbndaq-artdaq-core_Obj_SBND
		        sbndcode_CRT
		        sbndcode_CRTUtils
                    sbnobj::SBND_CRT
                    sbndcode_CRTUtils
                    sbndcode_GeoWrappers
                    nurandom::RandomUtils_NuRandomService_service
                    art::Framework_Core
                    art::Framework_Principal
                    art::Framework_Services_Registry
                    art::Framework_Services_Optional_RandomNumberGenerator_service
                    art::Persistency_Common
                    art::Persistency_Provenance
                    art::Utilities
                    canvas::canvas
                    messagefacility::MF_MessageLogger
                    fhiclcpp::fhiclcpp
                    cetlib::cetlib
                    CLHEP::CLHEP
                    ROOT::Core
                    ROOT::Tree
)


cet_build_plugin(HitDumper art::module SOURCE HitDumper_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(MuonTrackFilter art::module SOURCE MuonTrackFilter_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(MuonTrackProducer art::module SOURCE MuonTrackProducer_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(OnlineEvd art::module SOURCE OnlineEvd_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(TrackWiggles art::module SOURCE TrackWiggles_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(PMTRingingTagging art::module SOURCE PMTRingingTagging_module.cc LIBRARIES ${Ringing_LIBRARIES})
cet_build_plugin(GrabExplicitEventID art::module SOURCE GrabExplicitEventID_module.cc LIBRARIES ${MODULE_LIBRARIES})


install_source()

