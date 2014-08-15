#
#art_make_library( LIBRARY_NAME NuAna
#                  SOURCE NuAna_module
#                  LIBRARIES
#                        NuAna
#                        Geometry
#                        Geometry_service
#                        Simulation
#                        Utilities
#                        TimeService_service
#                        Filters
#                        RawData
#                        SignalShapingServiceT1034_service
#                        NuReweightArt
#                        ${SIMULATIONBASE}
#                        ${ART_FRAMEWORK_CORE}
#                        ${ART_FRAMEWORK_PRINCIPAL}
#                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
#                        ${ART_PERSISTENCY_COMMON}
#                        ${ART_PERSISTENCY_PROVENANCE}
#                        ${ART_UTILITIES}
#                        ${MF_MESSAGELOGGER}
#                        ${MF_UTILITIES}
#                        ${FHICLCPP}
#                        ${CETLIB}
#                        ${CLHEP}
#                        ${ROOT_BASIC_LIB_LIST}
#                  #BASENAME_ONLY
#)

# simple_plugin()

add_subdirectory(alg)
include_directories ( $ENV{GENIE_INC}/GENIE )

art_make( BASENAME_ONLY
          LIBRARY_NAME Ntuples          
          MODULE_LIBRARIES 
            NuAnaAlg
            Geometry_service
            TimeService_service
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
            Geometry
            Simulation
            Utilities
            Filters
            RawData
            NuReweightArt
            ${SIMULATIONBASE}
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
            ${ART_PERSISTENCY_COMMON}
            ${ART_PERSISTENCY_PROVENANCE}
            ${ART_UTILITIES}
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
            ${CLHEP}
            ${ROOT_BASIC_LIB_LIST}
          )

install_headers()
# install_fhicl()
# install_source()

