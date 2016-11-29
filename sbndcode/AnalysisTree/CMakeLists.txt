
art_make(
          MODULE_LIBRARIES larcore_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                           larsim_MCCheater_BackTracker_service
                           lardata_Utilities
                           larevt_Filters
                           lardataobj_RawData
                           lardataobj_RecoBase
                           lardata_RecoObjects larreco_Deprecated  
                           lardataobj_RecoBase
                           lardata_RecoObjects larreco_Deprecated  
                           larpandora_LArPandoraInterface
                           nusimdata_SimulationBase
                           ${LARRECO_LIB}
                           ${LARDATA_LIB}
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           art_Persistency_Common canvas_Persistency_Common
                           art_Persistency_Provenance canvas_Persistency_Provenance
                           art_Utilities canvas_Utilities
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CETLIB}
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
        )

# install_headers()
install_fhicl()
install_source()

