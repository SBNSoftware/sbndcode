art_dictionary(DICTIONARY_LIBRARIES 
lardataobj::RawData)

cet_build_plugin( GateTimingAnalyzer art::module
    SOURCE GateTimingAnalyzer_module.cc	
    LIBRARIES
    art_root_io::TFileService_service
    ROOT::Tree
    sbnobj::SBND_Timing
)

install_headers()
install_fhicl()
install_source()