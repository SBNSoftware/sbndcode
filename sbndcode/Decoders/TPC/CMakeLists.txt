cet_build_plugin( SBNDTPCDecoder art::module
  SOURCE SBNDTPCDecoder_module.cc
  LIBRARIES
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND_NevisTPC
  lardataobj::RawData
  artdaq_core::artdaq-core_Utilities
  art::Utilities
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
  art::Framework_Core
  ROOT::Core
  ROOT::Tree
  sbndcode_ChannelMaps_TPC_TPCChannelMapService_service
)

install_headers()
install_fhicl()
install_source()
art_dictionary()

