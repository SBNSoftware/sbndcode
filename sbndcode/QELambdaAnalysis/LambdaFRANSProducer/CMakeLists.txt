set(
MODULE_LIBRARIES
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art::Utilities
    canvas::canvas
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    art_root_io::tfile_support
    art_root_io::TFileService_service

    larcore::Geometry_Geometry_service
    larcorealg::Geometry 
    larevt::SpaceCharge
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service

    nusimdata::SimulationBase
    lardataobj::RecoBase
    larsim::MergeSimSources
    larsim::Utils
    larsim::MCCheater_ParticleInventoryService_service
    larsim::MCCheater_BackTrackerService_service

    ${ROOT_BASIC_LIB_LIST}
    ROOT::Tree
    ROOT::Core
    ROOT::TMVA

    sbndcode_FRANS
    sbndcode_TPCLines
    sbndcode_SObjects
    sbndcode_SEventHandle

    sbnobj::SBND_qelambda

    MichelVariables
    MichelExtractor
)

cet_build_plugin(FRANSProducer art::Module SOURCE FRANSProducer_module.cc LIBRARIES ${MODULE_LIBRARIES} )


install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)