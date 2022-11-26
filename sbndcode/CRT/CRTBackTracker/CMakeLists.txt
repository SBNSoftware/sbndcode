art_make_library(
                 LIBRARY_NAME sbndcode_CRT_CRTBackTracker
                 SOURCE       CRTBackTrackerAlg.cc
                 LIBRARIES
                        sbnobj_SBND_CRT
                        sbndcode_GeoWrappers
)

install_fhicl()
