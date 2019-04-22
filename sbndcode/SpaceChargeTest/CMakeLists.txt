simple_plugin(SpaceChargeTest "module"
                           larcore_Geometry_Geometry_service
                           lardata_Utilities
                           larevt_Filters
                           lardataobj_RawData
                           lardataobj_RecoBase
                           lardata_RecoObjects  
                           larpandora_LArPandoraInterface
                           nusimdata_SimulationBase
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                           art_Persistency_Common canvas
                           art_Persistency_Provenance canvas
                           art_Utilities canvas
                           ${MF_MESSAGELOGGER}
                           
                           ${FHICLCPP}
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
                           sbndcode_RecoUtils
			   sbndcode_SpaceCharge
        )

install_headers()
install_fhicl()
install_source()
