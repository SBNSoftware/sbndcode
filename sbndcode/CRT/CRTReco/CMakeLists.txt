art_make_library(
                 LIBRARY_NAME sbndcode_CRT_CRTReco
                 SOURCE       CRTClusterCharacterisationAlg.cc
                 LIBRARIES
                        sbnobj_SBND_CRT
                        sbndcode_GeoWrappers
)

simple_plugin(
              CRTStripHitProducer module
              larcorealg_Geometry
              sbnobj_SBND_CRT
              sbndcode_GeoWrappers
)

simple_plugin(
              CRTClusterProducer module
              larcorealg_Geometry
              sbnobj_SBND_CRT
              sbndcode_GeoWrappers
)

simple_plugin(
              CRTSpacePointProducer module
              larcorealg_Geometry
              sbnobj_SBND_CRT
              sbndcode_GeoWrappers
              sbndcode_CRT_CRTReco
)


install_fhicl()
