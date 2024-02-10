art_make_library(
                 SOURCE       TPCEventDisplayAlg.cc
                 LIBRARIES
                        larsim::Utils
                        ROOT::Graf3d
                        sbndcode_GeoWrappers
)

simple_plugin(
              TPCEventDisplay module
              sbndcode_TPCEventDisplay
)

install_fhicl()
