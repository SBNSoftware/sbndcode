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
  nusimdata_SimulationBase

  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}

  art_Persistency_Common
  canvas art_Persistency_Provenance
  canvas art_Utilities canvas
  ${MF_MESSAGELOGGER}

  ${FHICLCPP}
  ${ROOT_GEOM}
  ${ROOT_XMLIO}
  ${ROOT_GDML}
  ${ROOT_BASIC_LIB_LIST}

  sbndcode_RecoUtils
)

#install_headers()
install_fhicl()
install_source()
