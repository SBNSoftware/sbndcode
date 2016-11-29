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
#                        nutools_NuReweight_art
#                        nusimdata_SimulationBase
#                        ${ART_FRAMEWORK_CORE}
#                        ${ART_FRAMEWORK_PRINCIPAL}
#                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
#                        art_Persistency_Common canvas_Persistency_Common
#                        art_Persistency_Provenance canvas_Persistency_Provenance
#                        art_Utilities canvas_Utilities
#                        ${MF_MESSAGELOGGER}
#                        ${MF_UTILITIES}
#                        ${FHICLCPP}
#                        ${CETLIB}
#                        ${CLHEP}
#                        ${ROOT_BASIC_LIB_LIST}
#)

# simple_plugin()

add_subdirectory(alg)
include_directories ( $ENV{GENIE_INC}/GENIE )

art_make(
          LIBRARY_NAME Ntuples          
          MODULE_LIBRARIES 
            NuAnaAlg
            larcore_Geometry_Geometry_service
            # TimeService_service
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
            larcore_Geometry
            larsim_Simulation lardataobj_Simulation
            lardata_Utilities
            larevt_Filters
            lardataobj_RawData
            nutools_NuReweight_art
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
            ${ROOT_BASIC_LIB_LIST}
          )

install_headers()
install_fhicl()
# install_source()

