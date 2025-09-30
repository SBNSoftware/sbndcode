art_make_library(
                 SOURCE       CRTClusterCharacterisationAlg.cc
                 LIBRARIES
                 sbnobj::SBND_CRT
                 sbndcode::CRTGeoService
)

simple_plugin(
              CRTStripHitProducer module
              artdaq_core::artdaq-core_Data
              lardata::Utilities
              sbnobj::SBND_CRT
              sbnobj::SBND_Timing
              sbndcode::CRTGeoService
)

simple_plugin(
              CRTClusterProducer module
              lardata::Utilities
              sbnobj::SBND_CRT
              sbndcode::CRTGeoService
)

simple_plugin(
              CRTSpacePointProducer module
              lardata::Utilities
              sbnobj::SBND_CRT
              sbndcode_CRT_CRTReco
)

simple_plugin(
              CRTTrackProducer module
              lardata::Utilities
              sbnobj::SBND_CRT
              sbndcode_CRTUtils
              sbndcode::CRTGeoService
              Eigen3::Eigen
)

install_fhicl()
