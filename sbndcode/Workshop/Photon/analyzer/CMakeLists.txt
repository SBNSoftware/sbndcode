art_make(
MODULE_LIBRARIES
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_ROOT_IO_TFILE_SUPPORT}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas
  nusimdata_SimulationBase
  cetlib cetlib_except
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${ROOT_BASIC_LIB_LIST}

  lardataobj_Simulation

  sbndcode_OpDetSim
)
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
