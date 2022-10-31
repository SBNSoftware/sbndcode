art_make(
          MODULE_LIBRARIES
                           larpandora_LArPandoraInterface
                           lardata_RecoObjects  
                           lardataobj_RecoBase
                           lardataobj_AnalysisBase
                           larreco_RecoAlg

                           art_Persistency_Common canvas
                           art_Persistency_Provenance canvas
                           art_Utilities canvas
                           cetlib cetlib_except
                           ${LARRECO_LIB}
                           ${LARDATA_LIB}
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_BASIC}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CLHEP}
                           sbndcode_RecoUtils
                           sbndcode_OpDetSim

                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
                           ${ROOT_SPECTRUM}
                           ${ROOT_ROOFIT}
                           ${ROOT_ROOFITCORE}
        )
install_fhicl()
install_source()
