art_make(
MODULE_LIBRARIES
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  art_root_io::tfile_support
  art_root_io::TFileService_service
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  nusimdata::SimulationBase
  cetlib::cetlib cetlib_except
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ${ROOT_BASIC_LIB_LIST}

  lardataobj_Simulation

  sbndcode_OpDetSim
)
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
