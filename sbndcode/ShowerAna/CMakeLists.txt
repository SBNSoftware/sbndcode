set(
  MODULE_LIBRARIES
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  larsim::Simulation
  lardataobj::Simulation
  larsim::MCCheater_BackTrackerService_service
  lardata::Utilities
  larevt::Filters
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  lardataobj::MCBase
  larreco::RecoAlg
  lardata::RecoObjects
  larpandora::LArPandoraInterface
  nusimdata::SimulationBase
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Geom
  ROOT::XMLIO
  ROOT::Gdml
  ROOT::Core
  sbndcode_RecoUtils
)

cet_build_plugin(PFPSliceValidation art::module SOURCE PFPSliceValidation_module.cc LIBRARIES ${MODULE_LIBRARIES})

#install_headers()
install_fhicl()
install_source()
