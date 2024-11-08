cet_build_plugin( EventAna art::module
    SOURCE EventAna_module.cc
    LIBRARIES
    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND_NevisTPC
	sbnobj::SBND_Timing

    lardataobj::RecoBase
  
    art_root_io::TFileService_service
    art_root_io::tfile_support
	artdaq_core::artdaq-core_Data
    artdaq_core::artdaq-core_Utilities
    art::Utilities
    art::Framework_Core

    ROOT::Core
    ROOT::Tree
    ROOT::Spectrum

    fhiclcpp::fhiclcpp
    messagefacility::MF_MessageLogger
    )

install_fhicl()
