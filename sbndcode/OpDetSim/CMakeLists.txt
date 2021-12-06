
simple_plugin( sbndPDMapAlg tool
                    larcore_Geometry_Geometry_service
                    lardataobj_Simulation
                    lardata_Utilities
                    lardataobj_RawData
                    # lardataobj_RecoBase
                    lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
                    larpandora_LArPandoraInterface
                    sbndcode_Utilities_SignalShapingServiceSBND_service
                    nurandom::RandomUtils_NuRandomService_service
                    art::Framework_Services_Optional_RandomNumberGenerator_service
                    canvas::canvas
                    cetlib_except::cetlib_except
                    messagefacility::MF_MessageLogger
                    
                    fhiclcpp::fhiclcpp
                    cetlib::cetlib cetlib_except
                    CLHEP::CLHEP
                    ${ROOT_BASIC_LIB_LIST}
                    ROOT::Core
)


art_make( EXCLUDE sbndPDMapAlg_tool.cc
  LIB_LIBRARIES
                    sbndcode_OpDetSim_sbndPDMapAlg_tool
                    larcore_Geometry_Geometry_service
                    lardataobj_Simulation
                    lardata_Utilities
                    lardataobj_RawData
                    # lardataobj_RecoBase
                    lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
                    larpandora_LArPandoraInterface
                    sbndcode_Utilities_SignalShapingServiceSBND_service
                    nurandom::RandomUtils_NuRandomService_service
                    art::Framework_Services_Optional_RandomNumberGenerator_service
                    canvas::canvas
                    cetlib_except::cetlib_except
                    messagefacility::MF_MessageLogger
                    
                    fhiclcpp::fhiclcpp
                    cetlib::cetlib cetlib_except
                    CLHEP::CLHEP
                    ${ROOT_BASIC_LIB_LIST}
                    ROOT::Core

  MODULE_LIBRARIES
                    sbndcode_OpDetSim
                    pthread
                    larcore_Geometry_Geometry_service
                    lardataobj_Simulation
                    lardata_Utilities
                    lardataobj_RawData
                    lardataobj_RecoBase
                    lardataobj_AnalysisBase
                    lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
                    sbndcode_Utilities_SignalShapingServiceSBND_service
                    nurandom::RandomUtils_NuRandomService_service
                    art::Framework_Core
                    art::Framework_Principal
                    
                    art::Framework_Services_Optional_RandomNumberGenerator_service
                    canvas::canvas
                    messagefacility::MF_MessageLogger
                    
                    fhiclcpp::fhiclcpp
                    cetlib::cetlib cetlib_except
                    CLHEP::CLHEP
                    ${ROOT_BASIC_LIB_LIST}
)
install_headers()
install_fhicl()
install_source()
cet_enable_asserts()
add_subdirectory(FlashFinder)
add_subdirectory(PMTAlg)

# install sbnd_pds_mapping.json with mapping of the photon detectors
install_fw(LIST sbnd_pds_mapping.json)
