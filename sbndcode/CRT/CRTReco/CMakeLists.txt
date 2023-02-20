art_make_library(
                 LIBRARY_NAME sbndcode_CRT_CRTReco
                 SOURCE       CRTClusterCharacterisationAlg.cc
                 LIBRARIES
                        sbnobj_SBND_CRT
                        sbndcode_GeoWrappers
)

simple_plugin(
              CRTStripHitProducer module
              larcorealg::Geometry
              sbnobj::SBND_CRT
              sbndcode_GeoWrappers
)

simple_plugin(
              CRTClusterProducer module
              larcorealg::Geometry
              sbnobj::SBND_CRT
              sbndcode_GeoWrappers
)

simple_plugin(
              CRTSpacePointProducer module
              larcorealg::Geometry
              sbnobj::SBND_CRT
              sbndcode_GeoWrappers
              sbndcode_CRT_CRTReco
)

simple_plugin(
              CRTTrackProducer module
              larcorealg::Geometry
              lardata::Utilities
              sbnobj::SBND_CRT
              sbndcode_CRTUtils
              Eigen3::Eigen
)

install_fhicl()
