cet_build_plugin( TimingCAEN art::module
    SOURCE TimingCAEN_module.cc 
    LIBRARIES 
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
        sbndcode::ChannelMaps_PMT_PMTChannelMapService_service
        sbnobj::SBND_Timing
        sbnobj::SBND_CRT
        lardata::Utilities
        canvas::canvas
        lardataobj::RecoBase
        lardataobj::RawData
        lardataobj::headers
        art_root_io::TFileService_service
        art_root_io::tfile_support
        artdaq_core::artdaq-core_Data
        ROOT::Core
        ROOT::Tree
)

install_fhicl()
install_source()
install_headers()
