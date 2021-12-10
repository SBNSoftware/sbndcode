art_make(
  MODULE_LIBRARIES
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larsim_Simulation
  lardataobj_Simulation
  larsim_MCCheater_BackTrackerService_service
  lardata_Utilities
  larevt_Filters
  lardataobj_RawData
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  lardataobj_MCBase
  larreco_RecoAlg
  lardata_RecoObjects
  larpandora_LArPandoraInterface
  nusimdata::SimulationBase

  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  

  art::Persistency_Common
  canvas art::Persistency_Provenance
  canvas art::Utilities canvas
  messagefacility::MF_MessageLogger

  fhiclcpp::fhiclcpp
  ROOT::Geom
  ROOT::XMLIO
  ROOT::Gdml
  ${ROOT_BASIC_LIB_LIST}

  sbndcode_RecoUtils
)

#install_headers()
install_fhicl()
install_source()
