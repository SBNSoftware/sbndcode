art_make_library(
                 LIBRARY_NAME sbndcode_CRT_CRTTPCMatching
                 SOURCE       CRTTrackMatchAlg.cc
                              CRTSpacePointMatchAlg.cc
                 LIBRARIES
                        sbnobj_SBND_CRT
                        sbndcode_GeoWrappers
)

simple_plugin(
              CRTTrackMatching module
              sbnobj::SBND_CRT
              sbndcode_GeoWrappers
              sbndcode_CRT_CRTTPCMatching
)

simple_plugin(
              CRTSpacePointMatching module
              sbnobj::SBND_CRT
              sbndcode_GeoWrappers
              sbndcode_CRT_CRTTPCMatching
)

install_fhicl()
