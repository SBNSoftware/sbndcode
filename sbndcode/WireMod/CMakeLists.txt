find_library(sbncode::sbncode_WireMod_Utility_LIBRARY NAMES sbncode::sbncode_WireMod_Utility)

find_library(sbncode_CAFMaker_LIBRARY NAMES sbncode_CAFMaker)

if(NOT DEFINED sbncode_CAFMaker_LIBRARY)
    message(FATAL_ERROR "Could not find CAFMaker")
endif()

message(${sbncode_CAFMaker_LIBRARY})

find_package(sbnobj)

cet_build_plugin(WireModifier art::EDProducer
    SOURCE
      WireModifier_module.cc
      bspline.cpp
    LIBRARIES PRIVATE
    sbncode_CAFMaker
    # ${sbncode_CAFMaker_LIBRARY}
    sbncode::WireMod_Utility
    sbnobj::Common_SBNEventWeight
    larsim::MCCheater_BackTrackerService_service
    sbnanaobj::StandardRecord
)


install_headers()
install_fhicl()
install_source()
