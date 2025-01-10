add_subdirectory(PDSAnalyzer)

set(
MODULE_LIBRARIES
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities
  nusimdata::SimulationBase
  cetlib::cetlib
  cetlib_except::cetlib_except
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Core
  ROOT::Tree
  lardataobj::Simulation
  sbndcode_OpDetSim
)
cet_build_plugin(OpDetAnalyzer art::module SOURCE OpDetAnalyzer_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
