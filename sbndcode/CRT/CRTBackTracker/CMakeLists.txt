art_make_library(
                 SOURCE       CRTBackTrackerAlg.cc
                 LIBRARIES
                        sbnobj::SBND_CRT
                        sbndcode::TPCGeoAlg
                        sbndcode_Geometry_GeometryWrappers_CRTGeoService_service
                        larsim::MCCheater_ParticleInventoryService_service
)

install_fhicl()
