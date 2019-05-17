simple_plugin(SPhaseChannelNoiseService   "service"                                                           
                
		larcorealg_Geometry
		sbndcode_Utilities_SignalShapingServiceSBND_service
		${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
		${ART_ROOT_IO_TFILESERVICE_SERVICE}
		nutools_RandomUtils_NuRandomService_service
		${ART_FRAMEWORK_CORE}
		art_Utilities canvas
		cetlib cetlib_except
		${CLHEP}
 		${ROOT_BASIC_LIB_LIST}                                                                                                                   
)      


install_fhicl()
