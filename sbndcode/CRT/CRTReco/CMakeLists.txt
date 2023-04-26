art_make_library(
                 LIBRARY_NAME sbndcode_CRT_CRTReco
                 SOURCE       CRTClusterCharacterisationAlg.cc
                 LIBRARIES
                        sbnobj_SBND_CRT
                        sbndcode_GeoWrappers
)

simple_plugin(
              CRTStripHitProducer module
              sbnobj::SBND_CRT
              sbndcode_GeoWrappers
)

simple_plugin(
              CRTClusterProducer module
              sbnobj::SBND_CRT
              sbndcode_GeoWrappers
)

simple_plugin(
              CRTSpacePointProducer module
              sbnobj::SBND_CRT
              sbndcode_GeoWrappers
              sbndcode_CRT_CRTReco
)

simple_plugin(
              CRTTrackProducer module
              sbnobj::SBND_CRT
              sbndcode_CRTUtils
              sbndcode_GeoWrappers
              Eigen3::Eigen
)

install_fhicl()
