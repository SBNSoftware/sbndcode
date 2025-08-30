set (
  MODULE_LIBRARIES
                    lardata::Utilities
                    art_root_io::tfile_support
                    art_root_io::TFileService_service
                    artdaq_core::artdaq-core_Data
                    ROOT::Core
                    ROOT::Tree
                    sbnobj::SBND_Timing
                    sbnobj::SBND_CRT
)
cet_build_plugin(FrameShift art::module SOURCE FrameShift_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_fhicl()
install_source()
