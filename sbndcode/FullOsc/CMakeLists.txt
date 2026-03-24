set(
      MODULE_LIBRARIES
        PUBLIC
          art::Framework_Core
          art::Framework_Services_Registry
          art_root_io::TFileService_service
          art_root_io::TFileService_service
          art_root_io::tfile_support
          art::Framework_Principal
          art::Persistency_Common
          art::Utilities canvas::canvas
          cetlib::cetlib cetlib_except::cetlib_except
          CLHEP::CLHEP
          nugen::EventGeneratorBase_GENIE
          nusimdata::SimulationBase
          lardata::Utilities
          larcoreobj::SummaryData
          fhiclcpp::fhiclcpp
          messagefacility::MF_MessageLogger
)

cet_build_plugin(GENIEOscPileup art::module SOURCE GENIEOscPileup_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(GENIEFirst art::module SOURCE GENIEFirst_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()