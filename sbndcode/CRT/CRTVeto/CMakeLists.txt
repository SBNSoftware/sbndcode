art_make_library(
                 LIBRARIES
                        sbnobj::SBND_CRT
                        sbndcode_GeoWrappers
)

simple_plugin(
              CRTVetoProducer module
              artdaq_core::artdaq-core_Data
              sbnobj::SBND_Timing
              sbnobj::SBND_CRT
              sbndcode_CRTUtils
	      lardata::Utilities
	      art_root_io::TFileService_service
)

install_fhicl()
