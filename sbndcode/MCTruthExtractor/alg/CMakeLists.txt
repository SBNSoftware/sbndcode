include_directories ( $ENV{GENIE_INC}/GENIE )



art_make(
          LIBRARY_NAME NuAnaAlg
          LIB_LIBRARIES 
            larcore_Geometry_Geometry_service
            larcore_Geometry
            nutools_NuReweight_art
            nutools_NuReweight
            nusimdata_SimulationBase
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
            art_Persistency_Common canvas_Persistency_Common
            art_Persistency_Provenance canvas_Persistency_Provenance
            art_Utilities canvas_Utilities
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
            ${CLHEP}
            ${Boost_SYSTEM_LIBRARY}
            ${ROOT_BASIC_LIB_LIST}
          )

install_headers()
# install_fhicl()
# install_source()

