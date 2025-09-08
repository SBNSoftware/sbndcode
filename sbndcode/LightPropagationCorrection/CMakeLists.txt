
set (
  MODULE_LIBRARIES
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::tfile_support
    art_root_io::TFileService_service
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities

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
    nusimdata::SimulationBase
    nug4::ParticleNavigation

    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
    sbndaq_artdaq_core::sbndaq-artdaq-core_Obj_SBND
    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common

    sbndcode_OpDetReco_OpFlash_FlashFinder

    sbnobj::Common_Reco
    sbnobj::SBND_CRT
    sbnobj::SBND_Timing
    sbnobj::SBND_OpFlashTiming
    sbndcode_OpDetSim
)

cet_build_plugin(LightPropagationCorrection art::module SOURCE LightPropagationCorrection_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(LightPropagationCorrectionAna art::module SOURCE LightPropagationCorrectionAna_module.cc LIBRARIES ${MODULE_LIBRARIES})

build_dictionary(DICTIONARY_LIBRARIES
  lardataobj::RecoBase
  canvas::canvas
)
add_subdirectory(job)
install_fhicl()
install_source()
