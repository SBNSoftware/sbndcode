set(
	 MODULE_LIBRARIES
	 ROOT::Tree
	 sbnobj_SBND_CRT
	 sbndcode_GeoWrappers
)

cet_build_plugin(CRTSharpsAnalysis art::Module SOURCE CRTSharpsAnalysis_module.cc LIBRARIES ${MODULE_LIBRARIES})
