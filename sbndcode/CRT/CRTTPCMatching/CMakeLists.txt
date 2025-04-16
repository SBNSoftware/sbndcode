art_make_library(
                 SOURCE       CRTTrackMatchAlg.cc
                              CRTSpacePointMatchAlg.cc
                 LIBRARIES
                        sbnobj::SBND_CRT
                        sbndcode_GeoWrappers
)

simple_plugin(
              CRTTrackMatching module
              sbnobj::SBND_CRT
              sbndcode_GeoWrappers
              sbndcode_CRT_CRTTPCMatching
              sbndcode_CRT_CRTEventDisplay
)

simple_plugin(
              CRTSpacePointMatching module
              sbnobj::SBND_CRT
              sbndcode_GeoWrappers
              sbndcode_CRT_CRTTPCMatching
              sbndcode_CRT_CRTEventDisplay
)

install_fhicl()
