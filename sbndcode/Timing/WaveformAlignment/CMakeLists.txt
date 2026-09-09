set (
  MODULE_LIBRARIES
                    lardata::Utilities
                    lardataobj::RecoBase
                    art_root_io::tfile_support
                    art_root_io::TFileService_service
                    canvas::canvas
                    ROOT::Gpad
                    ROOT::Core
                    ROOT::Tree
                    sbnobj::SBND_Timing
                    sbnobj::SBND_CRT
)
cet_build_plugin(WaveformAlignment art::module SOURCE WaveformAlignment_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_fhicl()
install_source()
