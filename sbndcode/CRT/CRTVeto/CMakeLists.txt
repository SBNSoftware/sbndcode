

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
              sbndcode_GeoWrappers
)

install_fhicl()
