set(
  MODULE_LIBRARIES
                    sbndcode_OpDetSim
                                larcorealg::Geometry
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
                    art::Framework_Services_Registry
                    art_root_io::tfile_support
                    art_root_io::TFileService_service
                                artdaq_core::artdaq-core_Data
                    artdaq_core::artdaq-core_Utilities
                    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
                    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
                    sbndcode_CRTUtils
                    sbnobj::Common_CRT
                    sbnobj::SBND_CRT
                    sbnobj::SBND_Trigger
                    sbndcode_CRTUtils
                    sbndcode_GeoWrappers
                    art::Framework_Services_Optional_RandomNumberGenerator_service
                    canvas::canvas
                    messagefacility::MF_MessageLogger
                    fhiclcpp::fhiclcpp
                    cetlib::cetlib
                    CLHEP::CLHEP
                    ROOT::Core
                    ROOT::Tree
)

add_subdirectory(CRT)
add_subdirectory(PMT)

cet_build_plugin(ArtdaqFragmentProducer art::module SOURCE ArtdaqFragmentProducer_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(MetricProducer art::module SOURCE MetricProducer_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()
