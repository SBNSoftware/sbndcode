cet_enable_asserts()

art_make(
  BASENAME_ONLY
  LIB_LIBRARIES
                    larcore_Geometry_Geometry_service
                    lardataobj_Simulation
                    lardata_Utilities
                    lardataobj_RawData
                    # lardataobj_RecoBase
                    lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
                    larpandora_LArPandoraInterface
                    sbndcode_Utilities_SignalShapingServiceSBND_service
                    nurandom_RandomUtils_NuRandomService_service
                    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                    canvas
                    cetlib_except
                    ${MF_MESSAGELOGGER}
                    ${MF_UTILITIES}
                    ${FHICLCPP}
                    cetlib cetlib_except
                    ${CLHEP}
                    ${ROOT_BASIC_LIB_LIST}
                    ${ROOT_CORE}

  MODULE_LIBRARIES
                    sbndcode_OpDetSim
                    pthread
                    larcore_Geometry_Geometry_service
                    lardataobj_Simulation
                    lardata_Utilities
                    lardataobj_RawData
                    lardataobj_RecoBase
                    lardataobj_AnalysisBase
                    lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
                    sbndcode_Utilities_SignalShapingServiceSBND_service
                    nurandom_RandomUtils_NuRandomService_service
                    ${ART_FRAMEWORK_CORE}
                    ${ART_FRAMEWORK_PRINCIPAL}
                    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                    canvas
                    ${MF_MESSAGELOGGER}
                    ${MF_UTILITIES}
                    ${FHICLCPP}
                    cetlib cetlib_except
                    ${CLHEP}
                    ${ROOT_BASIC_LIB_LIST}
)
install_headers()
install_fhicl()
install_source()
# install sbnd_pds_mapping.json with mapping of the photon detectors
install_fw(LIST sbnd_pds_mapping.json)

# # install .dat files with test bench data
# file(GLOB dat_files *.dat)
# install_fw(LIST ${dat_files})

# simple_plugin(SimPMTSBND "module"
#                            larcore_Geometry_Geometry_service
#                            lardataobj_Simulation
#                            lardata_Utilities
#                            lardataobj_RawData
#                            lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
#                            sbndcode_Utilities_SignalShapingServiceSBND_service
#                            nurandom_RandomUtils_NuRandomService_service
# 			   ${ART_FRAMEWORK_CORE}
#                            ${ART_FRAMEWORK_PRINCIPAL}
#                            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
#                            canvas
#                            ${MF_MESSAGELOGGER}
#                            ${MF_UTILITIES}
#                            ${FHICLCPP}
#                            cetlib cetlib_except
#                            ${CLHEP}
#                            ${ROOT_BASIC_LIB_LIST}
# )
#
# simple_plugin(SimArapucaSBND "module"
#                            larcore_Geometry_Geometry_service
#                            lardataobj_Simulation
#                            lardata_Utilities
#                            lardataobj_RawData
#                            lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
#                            sbndcode_Utilities_SignalShapingServiceSBND_service
#                            nurandom_RandomUtils_NuRandomService_service
#                            ${ART_FRAMEWORK_CORE}
#                            ${ART_FRAMEWORK_PRINCIPAL}
#                            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
#                            canvas
#                            ${MF_MESSAGELOGGER}
#                            ${MF_UTILITIES}
#                            ${FHICLCPP}
#                            cetlib cetlib_except
#                            ${CLHEP}
#                            ${ROOT_BASIC_LIB_LIST}
# )
