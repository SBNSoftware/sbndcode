art_make_library(
                 SOURCE       CRTBackTrackerAlg.cc
                 LIBRARIES
                        sbnobj::SBND_CRT
                        sbndcode_GeoWrappers
                        larsim::MCCheater_ParticleInventoryService_service
)

install_fhicl()
