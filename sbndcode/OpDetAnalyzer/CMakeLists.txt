set(
MODULE_LIBRARIES
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  nusimdata::SimulationBase
  cetlib::cetlib
  cetlib_except::cetlib_except
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Core
  lardataobj::Simulation
  sbndcode_OpDetSim
) 
cet_build_plugin(OpDetAnalyzer art::module SOURCE OpDetAnalyzer_module.cc LIBRARIES ${art::Framework_Core})

install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
