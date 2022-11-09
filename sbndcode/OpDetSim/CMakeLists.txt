
cet_build_plugin( sbndPDMapAlg art::tool
                  SOURCE
                    sbndPDMapAlg_tool.cc
                  LIBRARIES
                    larcore::Geometry_Geometry_service
                    lardataobj::Simulation
                    lardata::Utilities
                    lardataobj::RawData
                    # lardataobj::RecoBase
                    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
                    larpandora::LArPandoraInterface
                    sbndcode_Utilities_SignalShapingServiceSBND_service
                    nurandom::RandomUtils_NuRandomService_service
                    art::Framework_Services_Optional_RandomNumberGenerator_service
                    canvas::canvas
                    cetlib_except::cetlib_except
                    messagefacility::MF_MessageLogger
                    fhiclcpp::fhiclcpp
                    cetlib::cetlib
                    cetlib_except::cetlib_except
                    CLHEP::CLHEP
                    ROOT::Core
)

art_make_library( SOURCE
                    DigiArapucaSBNDAlg.cc
                    DigiPMTSBNDAlg.cc
                    opDetDigitizerSBND_module.cc
                    opDetDigitizerWorker.cc
                    opDetSBNDTriggerAlg.cc
                  LIBRARIES
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
                    ROOT::Core
)
set (
  MODULE_LIBRARIES
                    sbndcode_OpDetSim
                    larcore::Geometry_Geometry_service
                    lardataobj::Simulation
                    lardata::Utilities
                    lardataobj::RawData
                    lardataobj::RecoBase
                    lardataobj::AnalysisBase
                    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
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
                    ROOT::Core
)


cet_build_plugin(opHitFinderSBND art::module SOURCE opHitFinderSBND_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(opDetDigitizerSBND art::module SOURCE opDetDigitizerSBND_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(wvfAna art::module SOURCE wvfAna_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()
cet_enable_asserts()
add_subdirectory(FlashFinder)
add_subdirectory(PMTAlg)

# install sbnd_pds_mapping.json with mapping of the photon detectors
install_fw(LIST sbnd_pds_mapping.json)
