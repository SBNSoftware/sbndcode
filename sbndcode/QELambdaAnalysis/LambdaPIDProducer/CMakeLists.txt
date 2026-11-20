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

    sbnobj::SBND_qelambda
)

cet_build_plugin(LambdaPID art::Module SOURCE LambdaPID_module.cc LIBRARIES ${MODULE_LIBRARIES} )


install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)