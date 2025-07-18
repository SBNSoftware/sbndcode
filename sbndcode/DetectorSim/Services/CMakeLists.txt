cet_build_plugin(SBNDNoNoiseService   art::service
                 SOURCE SBNDNoNoiseService_service.cc
)

cet_build_plugin(SBNDThermalNoiseServiceInTime   art::service
              SOURCE SBNDThermalNoiseServiceInTime_service.cc
              LIBRARIES
		larcorealg::Geometry
		sbndcode_Utilities_SignalShapingServiceSBND_service
		art_root_io::tfile_support
		art_root_io::TFileService_service
		nurandom::RandomUtils_NuRandomService_service
		art::Framework_Core
		art::Utilities
		cetlib::cetlib
		CLHEP::CLHEP
 		ROOT::Core
)

cet_build_plugin(SBNDThermalNoiseServiceInFreq   art::service
               SOURCE SBNDThermalNoiseServiceInFreq_service.cc
               LIBRARIES
		larcorealg::Geometry
		sbndcode_Utilities_SignalShapingServiceSBND_service
		art_root_io::tfile_support ROOT::Core
		art_root_io::TFileService_service
		nurandom::RandomUtils_NuRandomService_service
		art::Framework_Core
		CLHEP::CLHEP
 		ROOT::Core
)

cet_build_plugin(SBNDuBooNEDataDrivenNoiseService   art::service
                SOURCE SBNDuBooNEDataDrivenNoiseService_service.cc
                LIBRARIES
		larcorealg::Geometry
		sbndcode_Utilities_SignalShapingServiceSBND_service
		art_root_io::tfile_support
		art_root_io::TFileService_service
		nurandom::RandomUtils_NuRandomService_service
		art::Framework_Core
		CLHEP::CLHEP
 		ROOT::Core
)

install_fhicl()
