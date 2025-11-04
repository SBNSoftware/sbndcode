art_make_library(
                 SOURCE       CRTEventDisplayAlg.cc
                 LIBRARIES
                        larsim::Utils
                        ROOT::Graf3d
                        sbnobj::SBND_CRT
                        sbndcode_CRT_CRTBackTracker
                        lardata::DetectorPropertiesService
)

simple_plugin(
              CRTEventDisplay module
              sbndcode_CRT_CRTEventDisplay
)

install_fhicl()
