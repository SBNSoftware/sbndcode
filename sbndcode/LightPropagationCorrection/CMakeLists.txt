
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
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
    larcore::Geometry_Geometry_service
    larcorealg::Geometry
    lardataobj::RecoBase
    lardataobj::MCBase
    lardataobj::Simulation


    sbndcode_OpDetReco_OpFlash_FlashFinder

    sbnobj::Common_Reco
    sbnobj::SBND_Timing
    sbndcode_OpDetSim
)

cet_build_plugin(LightPropagationCorrection art::module SOURCE LightPropagationCorrection_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(LightPropagationCorrectionAna art::module SOURCE LightPropagationCorrectionAna_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_fhicl()
install_source()
install_headers()
add_subdirectory(job)

