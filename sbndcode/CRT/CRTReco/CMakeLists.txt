art_make_library(
                 SOURCE       CRTClusterCharacterisationAlg.cc
                 LIBRARIES
                 sbnobj::SBND_CRT
                 sbndcode_Geometry_GeometryWrappers_CRTGeoService_service
)

cet_build_plugin(
              CRTStripHitProducer art::module
              LIBRARIES
              artdaq_core::artdaq-core_Data
              lardata::Utilities
              sbnobj::SBND_CRT
              sbnobj::SBND_Timing
              sbndcode_Geometry_GeometryWrappers_CRTGeoService_service
)

cet_build_plugin(
              CRTClusterProducer art::module
              LIBRARIES
              lardata::Utilities
              sbnobj::SBND_CRT
              sbndcode_Geometry_GeometryWrappers_CRTGeoService_service
)

cet_build_plugin(
              CRTSpacePointProducer art::module
              LIBRARIES
              lardata::Utilities
              sbnobj::SBND_CRT
              sbndcode_CRT_CRTReco
)

cet_build_plugin(
              CRTTrackProducer art::module
              LIBRARIES
              lardata::Utilities
              sbnobj::SBND_CRT
              sbndcode_CRTUtils
              sbndcode_Geometry_GeometryWrappers_CRTGeoService_service
              Eigen3::Eigen
)

cet_build_plugin(
              CRTBlobProducer art::module
              LIBRARIES
              lardata::Utilities
              sbnobj::SBND_CRT
)

install_fhicl()
