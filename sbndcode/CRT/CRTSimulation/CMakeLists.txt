art_make_library(
                 SOURCE       CRTDetSimAlg.cxx
                 LIBRARIES
                        ROOT::MathMore
                        fhiclcpp::fhiclcpp
                        art::Persistency_Provenance
                        canvas::canvas
                        sbnobj::SBND_CRT
                        sbndcode::CRTGeoService
)

simple_plugin(
              CRTDetSim module
              nurandom::RandomUtils_NuRandomService_service
              larcorealg::Geometry
              lardata::DetectorClocksService
              sbnobj::SBND_CRT
              sbndcode_CRT_CRTSimulation
)

install_headers()
install_source()
install_fhicl()
