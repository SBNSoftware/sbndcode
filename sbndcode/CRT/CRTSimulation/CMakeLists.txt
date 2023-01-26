art_make_library(
                 LIBRARY_NAME sbndcode_CRT_CRTSimulation
                 SOURCE       CRTDetSimAlg.cxx
                 LIBRARIES
                        ROOT::MathMore
                        fhiclcpp::fhiclcpp
                        art::Persistency_Provenance
                        canvas::canvas
                        sbnobj_SBND_CRT
                        sbndcode_GeoWrappers
)

simple_plugin(
              CRTDetSim module
              nurandom::RandomUtils_NuRandomService_service
              larcorealg::Geometry
              sbnobj::SBND_CRT
              sbndcode_CRT_CRTSimulation
)

install_headers()
install_source()
install_fhicl()
