cet_build_plugin( BeamAnalysis art::module
    SOURCE BeamAnalysis_module.cc 
    LIBRARIES 
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
        sbndaq_artdaq_core::sbndaq-artdaq-core_Obj_SBND
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
