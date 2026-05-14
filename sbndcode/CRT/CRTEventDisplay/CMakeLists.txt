art_make_library(
                 SOURCE       CRTEventDisplayAlg.cc
                 LIBRARIES
                        larsim::Utils
                        ROOT::Graf3d
                        sbnobj::SBND_CRT
                        sbndcode_CRT_CRTBackTracker
                        lardata::DetectorPropertiesService
)

cet_build_plugin(
              CRTEventDisplay art::module
              LIBRARIES sbndcode_CRT_CRTEventDisplay
)

cet_build_plugin(
              CRTChannelMappingEventDisplay art::module
              LIBRARIES sbndcode_CRT_CRTEventDisplay
)

install_fhicl()
