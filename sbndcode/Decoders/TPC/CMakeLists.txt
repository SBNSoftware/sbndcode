
# this is borrowed from sbndaq-artdaq-core to make the c7 compiler happy

add_definitions(-Wno-nested-anon-types)

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

cet_build_plugin( NTBDecoder art::module
  SOURCE NTBDecoder_module.cc
  LIBRARIES
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
  artdaq_core::artdaq-core_Utilities
  art::Utilities
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
  art::Framework_Core
  ROOT::Core
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()
art_dictionary()

