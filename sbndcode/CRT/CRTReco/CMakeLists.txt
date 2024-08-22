art_make_library(
                 SOURCE       CRTClusterCharacterisationAlg.cc
                 LIBRARIES
                        sbnobj::SBND_CRT
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
