install_headers()
install_fhicl()
install_source()

file(GLOB channel_map_file *ChannelMap*.txt)
install_fw( LIST ${channel_map_file} )

art_make(SERVICE_LIBRARIES 
                           art::Framework_Services_Registry
                           art::Framework_Principal
                           art::Framework_Core
                           art::Persistency_Provenance
                           messagefacility::MF_MessageLogger
                           ROOT::Core
         MODULE_LIBRARIES
                           art_root_io::TFileService_service
                           lardataobj::RawData
                           larcore::Geometry_Geometry_service
                           sbndcode_ChannelMaps_TPC_TPCChannelMapService_service
                           ROOT::Tree
)


