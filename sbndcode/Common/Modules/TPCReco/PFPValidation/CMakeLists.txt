art_make(MODULE_LIBRARIES
          sbndcode_Common_Modules_MCRecoUtils
          lardataobj_RecoBase
          larsim_Utils
	  lardataalg_DetectorInfo
          ${ART_FRAMEWORK_SERVICES_REGISTRY}
          ROOT::Core
          ROOT::Physics
          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
          ${MF_MESSAGELOGGER}
          ${ART_ROOT_IO_TFILE_SUPPORT}
        )

install_headers()
install_fhicl()
install_source()
