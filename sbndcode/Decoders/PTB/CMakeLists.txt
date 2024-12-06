
# this is borrowed from sbndaq-artdaq-core to make the c7 compiler happy
# remove if not needed by the trigger decoder
# add_definitions(-Wno-nested-anon-types)

cet_build_plugin( SBNDPTBDecoder art::module
  SOURCE SBNDPTBDecoder_module.cc
  LIBRARIES
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
  lardataobj::RawData
  artdaq_core::artdaq-core_Utilities
  art::Utilities
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
  art::Framework_Core
  ROOT::Core
)

cet_make_library(SOURCE
  SBNDPTBRawUtils.cxx
  LIBRARIES
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
  ROOT::Core
)

install_headers()
install_fhicl()
install_source()
art_dictionary()

