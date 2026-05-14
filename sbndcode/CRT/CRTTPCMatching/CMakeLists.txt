art_make_library(
                 SOURCE       CRTTrackMatchAlg.cc
                              CRTSpacePointMatchAlg.cc
                 LIBRARIES
                        sbnobj::SBND_CRT
                        sbndcode_Geometry_GeometryWrappers_CRTGeoService_service
                        lardata::DetectorPropertiesService
                        lardataobj::AnalysisBase
)

cet_build_plugin(
              CRTTrackMatching art::module
              LIBRARIES
              sbndcode_CRT_CRTTPCMatching
              sbndcode_CRT_CRTEventDisplay
)

cet_build_plugin(
              CRTSpacePointMatching art::module
              LIBRARIES
              sbndcode_CRT_CRTTPCMatching
              sbndcode_CRT_CRTEventDisplay
)

install_fhicl()
