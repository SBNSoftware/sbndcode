add_definitions(-Wno-nested-anon-types)

cet_build_plugin( FragmentScrambler art::module
  SOURCE ArtDAQFragmentBlender_module.cc
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
  sbnobj::Common_Analysis
)
find_package(gallery         REQUIRED)
find_package(lardataalg      REQUIRED)
find_package(lardataobj      REQUIRED)
find_package(larcorealg      REQUIRED)
find_package(larcoreobj      REQUIRED)
find_package(sbnanaobj       REQUIRED)
find_package(sbnobj          REQUIRED)
find_package(sbncode         REQUIRED)
find_package(sbndcode        REQUIRED)
find_package(fhiclcpp        REQUIRED)
find_package(messagefacility REQUIRED)
find_package(ROOT
  COMPONENTS Hist Tree RIO
  REQUIRED
  )


install_headers()
install_fhicl()
install_source()
art_dictionary()