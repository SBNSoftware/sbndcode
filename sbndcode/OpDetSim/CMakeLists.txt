
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
                    sbncode::OpDet_PDMapAlgSimple_tool
                    sbndcode_Utilities_SignalShapingServiceSBND_service
                    nurandom::RandomUtils_NuRandomService_service
                    art::Framework_Services_Optional_RandomNumberGenerator_service
                    canvas::canvas
                    messagefacility::MF_MessageLogger
                    fhiclcpp::fhiclcpp
                    cetlib::cetlib
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
                    larcore::Geometry_Geometry_service
                    lardataobj::Simulation
                    lardata::Utilities
                    lardataobj::RawData
                    # lardataobj::RecoBase
                    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
                    larpandora::LArPandoraInterface
                    sbndcode::Utilities_SignalShapingServiceSBND_service
                    nurandom::RandomUtils_NuRandomService_service
                    art::Framework_Services_Optional_RandomNumberGenerator_service
                    canvas::canvas
                    messagefacility::MF_MessageLogger
                    fhiclcpp::fhiclcpp
                    cetlib::cetlib
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
                    art::Framework_Services_Optional_RandomNumberGenerator_service
                    canvas::canvas
                    messagefacility::MF_MessageLogger
                    fhiclcpp::fhiclcpp
                    cetlib::cetlib
                    CLHEP::CLHEP
                    ROOT::Core
)

set (
  OTHER_LIBRARIES
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
                    cetlib::cetlib
                    CLHEP::CLHEP
                    ROOT::Core
                    larcorealg::Geometry
                    larcore::Geometry_Geometry_service
                    larsim::Simulation
                    lardataobj::Simulation
                    larsim::MCCheater_BackTrackerService_service
                    lardata::Utilities
                    larevt::Filters
                    lardataobj::RawData
                    lardataobj::RecoBase
                    lardataobj::AnalysisBase
                    lardataobj::MCBase
                    larreco::RecoAlg
                    lardata::RecoObjects
                    larpandora::LArPandoraInterface
                    nusimdata::SimulationBase
                    art::Framework_Core
                    art::Framework_Principal
                    art::Framework_Services_Registry
                    art::Persistency_Common
                    art::Persistency_Provenance
                    art::Utilities
                    messagefacility::MF_MessageLogger
                    fhiclcpp::fhiclcpp
                    ROOT::Geom
                    ROOT::XMLIO
                    ROOT::Gdml
                    ROOT::Core
                    sbndcode_RecoUtils
)


cet_build_plugin(opHitFinderSBND art::module SOURCE opHitFinderSBND_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(opDetDigitizerSBND art::module SOURCE opDetDigitizerSBND_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(wvfAna art::module SOURCE wvfAna_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(BeamRatesCalib art::module SOURCE BeamRatesCalib_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(PMTNoiseCounter art::module SOURCE PMTNoiseCounter_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(FlashWithNeutrino art::module SOURCE FlashWithNeutrinoInTPC.cc LIBRARIES ${OTHER_LIBRARIES})
cet_build_plugin(SummedWaveformProducer art::module SOURCE Zero_Surpressed_PMT_Summer_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()
cet_enable_asserts()
add_subdirectory(PMTAlg)
add_subdirectory(HDWvf)

# install sbnd_pds_mapping.json with mapping of the photon detectors
install_fw(LIST sbnd_pds_mapping.json)
