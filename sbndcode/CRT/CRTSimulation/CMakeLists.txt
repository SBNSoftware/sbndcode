art_make_library(
                 SOURCE       CRTDetSimAlg.cxx
                 LIBRARIES
                        ROOT::MathMore
                        fhiclcpp::fhiclcpp
                        art::Persistency_Provenance
                        canvas::canvas
                        sbnobj::SBND_CRT
                        sbndcode_Geometry_GeometryWrappers_CRTGeoService_service
)

cet_build_plugin(
              CRTDetSim art::module
              LIBRARIES
              nurandom::RandomUtils_NuRandomService_service
              larcorealg::Geometry
              lardata::DetectorClocksService
              sbnobj::SBND_CRT
              sbndcode_CRT_CRTSimulation
)

install_headers()
install_source()
install_fhicl()
