
art_make(
          MODULE_LIBRARIES larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation lardataobj_Simulation
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           lardata_Utilities
                           larevt_Filters
                           lardataobj_RawData
                           lardataobj_RecoBase
                           lardataobj_AnalysisBase
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
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
                           ${ROOT_SPECTRUM}
                           ${ROOT_ROOFIT}
                           ${ROOT_ROOFITCORE} 
                           sbndcode_RecoUtils
                           sbndcode_CRT
                           sbndcode_CRTProducts
                           sbndcode_CRTUtils
                           sbndcode_GeoWrappers
        )

install_headers()
install_fhicl()
install_source()

