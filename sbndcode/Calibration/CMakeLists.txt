add_subdirectory(tools)

art_make(MODULE_LIBRARIES
			larcore_Geometry_Geometry_service
                        lardata_Utilities
			lardata_ArtDataHelper
                        lardataobj_RawData
                        lardataobj_RecoBase
                        sbndcode_Utilities_SignalShapingServiceT1053_service
                        sbndcode_Utilities_SignalShapingServiceSBND_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			art_Utilities			
                        canvas
                        ${MF_MESSAGELOGGER}
                        
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}

TOOL_LIBRARIES
			larcore_Geometry_Geometry_service
			lardata_Utilities
                        lardata_ArtDataHelper
                        lardataobj_RawData
                        lardataobj_RecoBase
                        sbndcode_Utilities_SignalShapingServiceSBND_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			art_Utilities			
                        canvas
                        ${MF_MESSAGELOGGER}

                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()
