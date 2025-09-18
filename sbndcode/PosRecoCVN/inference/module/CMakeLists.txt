set(MODULE_LIBRARIES
    art::Framework_Principal
    art::Framework_Services_Registry
    art::Framework_Core
    art_root_io::tfile_support
    art_root_io::TFileService_service
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities
    canvas::canvas
    messagefacility::MF_MessageLogger

    ROOT::Tree
    ROOT::Core
    
    larsim::Utils
    larsim::MCCheater_BackTrackerService_service
    larsim::MCCheater_ParticleInventoryService_service
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
    larcore::Geometry_Geometry_service
    larcorealg::Geometry
    lardataobj::RawData
    lardataobj::RecoBase
    lardataobj::MCBase
    lardataobj::Simulation
    lardataobj::AnalysisBase
    nusimdata::SimulationBase
    nug4::ParticleNavigation

    sbnobj::Common_Reco
    sbndcode_OpDetSim
    sbndcode_PosRecoCVN_inference_tf

)

cet_build_plugin(PosRecoCVNProducer art::Module SOURCE PosRecoCVNProducer_module.cc LIBRARIES ${MODULE_LIBRARIES} )

# Build ROOT dictionary
art_dictionary(DICTIONARY_LIBRARIES
    canvas::canvas
    ROOT::Core
)

install_fhicl()
install_source()