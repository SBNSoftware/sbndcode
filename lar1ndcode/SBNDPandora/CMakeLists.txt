
# where should the scripts/..xml file be installed?  Perhaps in bin?

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

art_make( BASENAME_ONLY
          LIBRARY_NAME  SBNDPandora
          LIB_LIBRARIES larpandora_LArPandoraInterface
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation
                        ${LARCONTENT}
                        ${PANDORASDK}
                        ${PANDORAMONITORING}
                        ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${Boost_SYSTEM_LIBRARY}
          MODULE_LIBRARIES larpandora_LArPandoraInterface 
                           SBNDPandora
          )

install_headers()
install_fhicl()
install_source()

add_subdirectory(scripts)
