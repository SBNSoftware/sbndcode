add_subdirectory(Analyzer)

art_make( BASENAME_ONLY MODULE_LIBRARIES
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_RawData
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_PhotonBackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  larsim_Simulation
  nug4_ParticleNavigation
  lardataobj_Simulation
  nusimdata_SimulationBase	
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${ROOT_GEOM}
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  )

install_fhicl()
install_headers()
install_source()
