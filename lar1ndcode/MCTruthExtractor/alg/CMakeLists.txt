
art_make( BASENAME_ONLY
          LIBRARY_NAME NuAnaAlg
          LIB_LIBRARIES 
            Geometry_service
            Geometry
            # Simulation
            # Utilities
            # Filters
            # RawData
            NuReweightArt
            NuReweight
            ${SIMULATIONBASE}
            ${ART_FRAMEWORK_CORE}
            # ${ART_FRAMEWORK_PRINCIPAL}
            # ${ART_FRAMEWORK_SERVICES_REGISTRY}
            # ${ART_FRAMEWORK_SERVICES_OPTIONAL}
            # ${ART_PERSISTENCY_COMMON}
            # ${ART_PERSISTENCY_PROVENANCE}
            # ${ART_UTILITIES}
            # ${MF_MESSAGELOGGER}
            # ${MF_UTILITIES}
            # ${FHICLCPP}
            # ${CETLIB}
            # ${CLHEP}
            ${ROOT_BASIC_LIB_LIST}
          )

install_headers()
# install_fhicl()
# install_source()

