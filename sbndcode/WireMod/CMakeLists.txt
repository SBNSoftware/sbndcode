find_package(sbnobj)

cet_build_plugin(WireModifier art::EDProducer
    SOURCE
      WireModifier_module.cc
    LIBRARIES PRIVATE
    sbncode_CAFMaker
    sbncode::WireMod_Utility
    sbnobj::Common_SBNEventWeight
    larsim::MCCheater_BackTrackerService_service
    sbnanaobj::StandardRecord
)


install_headers()
install_fhicl()
install_source()
