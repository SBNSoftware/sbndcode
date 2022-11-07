art_make_library(
		 LIBRARY_NAME sbndcode_CRT_CRTEventDisplay
		 SOURCE	      CRTEventDisplayAlg.cc
		 LIBRARIES
			larsim_Utils
	 		sbnobj_SBND_CRT
	 		sbndcode_GeoWrappers
)

simple_plugin(
	      CRTEventDisplay module
	      sbndcode_CRT_CRTEventDisplay
)

install_fhicl()
