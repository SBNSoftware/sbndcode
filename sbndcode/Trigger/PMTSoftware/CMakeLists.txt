set(
  MODULE_LIBRARIES
                    sbndcode_OpDetSim
                    lardataobj::RawData
                    art::Framework_Core
                    art::Framework_Principal
                    art::Framework_Services_Registry
                    art_root_io::tfile_support
                    art_root_io::TFileService_service
                    artdaq_core::artdaq-core_Data
                    artdaq_core::artdaq-core_Utilities
                    sbnobj::SBND_Timing
                    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
                    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
                    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND 
                    sbndaq_artdaq_core::sbndaq-artdaq-core_Obj_SBND
                    canvas::canvas
                    messagefacility::MF_MessageLogger
                    fhiclcpp::fhiclcpp
                    cetlib::cetlib
                    CLHEP::CLHEP
                    ROOT::Core
                    ROOT::Tree
)

cet_build_plugin(MetricAnalyzer art::module SOURCE MetricAnalyzer_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(MetricFilter   art::module SOURCE MetricFilter_module.cc   LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(PMTMetricProducer art::module SOURCE PMTMetricProducer_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(PMTMCMetricProducer art::module SOURCE PMTMCMetricProducer_module.cc LIBRARIES ${MODULE_LIBRARIES})

add_subdirectory(fcls)

install_headers()
install_fhicl()
install_source()