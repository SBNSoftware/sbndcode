simple_plugin(SBNDNoNoiseService   "service"                                                           
                
)

simple_plugin(SBNDThermalNoiseServiceInTime   "service"                                                           
                
		larcorealg_Geometry
		sbndcode_Utilities_SignalShapingServiceSBND_service
		art_root_io::tfile_support ROOT::Core
		art_root_io::TFileService_service
		nurandom::RandomUtils_NuRandomService_service
		art::Framework_Core
		art::Utilities canvas
		cetlib::cetlib cetlib_except
		CLHEP::CLHEP
 		${ROOT_BASIC_LIB_LIST}                                                                                                                   
)

simple_plugin(SBNDThermalNoiseServiceInFreq   "service"               
                
		larcorealg_Geometry
		sbndcode_Utilities_SignalShapingServiceSBND_service
		art_root_io::tfile_support ROOT::Core
		art_root_io::TFileService_service
		nurandom::RandomUtils_NuRandomService_service
		art::Framework_Core
		art::Utilities canvas
		cetlib::cetlib cetlib_except
		CLHEP::CLHEP
 		${ROOT_BASIC_LIB_LIST}                                                                                                                   
)      

simple_plugin(SBNDuBooNEDataDrivenNoiseService   "service"                                                           
                
		larcorealg_Geometry
		sbndcode_Utilities_SignalShapingServiceSBND_service
		art_root_io::tfile_support ROOT::Core
		art_root_io::TFileService_service
		nurandom::RandomUtils_NuRandomService_service
		art::Framework_Core
		art::Utilities canvas
		cetlib::cetlib cetlib_except
		CLHEP::CLHEP
 		${ROOT_BASIC_LIB_LIST}                                                                                                                   
)      


install_fhicl()
