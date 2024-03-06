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
)


cet_build_plugin(HitDumper art::module SOURCE HitDumper_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(MuonTrackFilter art::module SOURCE MuonTrackFilter_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(MuonTrackProducer art::module SOURCE MuonTrackProducer_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_source()

