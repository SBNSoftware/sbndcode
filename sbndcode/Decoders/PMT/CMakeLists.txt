art_make(
	MODULE_LIBRARIES 
	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
	sbnobj::SBND_Timing

    lardataobj::RecoBase
  
    art_root_io::TFileService_service
    art_root_io::tfile_support
	artdaq_core::artdaq-core_Data

    ROOT::Core
    ROOT::Tree

)

install_fhicl()
