cet_build_plugin( SBNDPMTDecoder art::module
    SOURCE SBNDPMTDecoder_module.cc	
    LIBRARIES 
	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
	sbnobj::SBND_Timing
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

art_dictionary(DICTIONARY_LIBRARIES 
                lardataobj::RawData)

install_fhicl()
install_headers()
install_source()
