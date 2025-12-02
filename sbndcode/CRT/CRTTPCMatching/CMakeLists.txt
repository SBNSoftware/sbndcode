art_make_library(
                 SOURCE       CRTTrackMatchAlg.cc
                              CRTSpacePointMatchAlg.cc
                 LIBRARIES
                        sbnobj::SBND_CRT
                        sbndcode_Geometry_GeometryWrappers_CRTGeoService_service
                        lardata::DetectorPropertiesService
                        lardataobj::AnalysisBase
)

simple_plugin(
              CRTTrackMatching module
              sbndcode_CRT_CRTTPCMatching
              sbndcode_CRT_CRTEventDisplay
)

simple_plugin(
              CRTSpacePointMatching module
              sbndcode_CRT_CRTTPCMatching
              sbndcode_CRT_CRTEventDisplay
)

install_fhicl()
