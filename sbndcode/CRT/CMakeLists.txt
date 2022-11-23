add_subdirectory(CRTAna)
add_subdirectory(CRTEventDisplay)
add_subdirectory(CRTReco)
add_subdirectory(CRTSimulation)
add_subdirectory(CRTUtils)

art_make(
	 EXCLUDE
		CRTChannelMapAlg.cxx
		CRTGeometryHelper_service.cc
	 LIBRARY_NAME
		sbndcode_CRTData
	 LIB_LIBRARIES
		larcorealg_Geometry
		larcore_Geometry_Geometry_service
		lardataobj_RecoBase
		art::Framework_Core
		art::Framework_IO_Sources
		art::Framework_Principal
		canvas::canvas
		art::Framework_Services_Registry
		messagefacility::MF_MessageLogger
		
		fhiclcpp::fhiclcpp
		${ROOT_BASIC_LIB_LIST}
		ROOT::Geom
		cetlib::cetlib
)

art_make_library(
		 LIBRARY_NAME sbndcode_CRT
    		 SOURCE       CRTChannelMapAlg.cxx
		 LIBRARIES
			larcorealg_Geometry
        		sbndcode_CRTData
        		cetlib_except::cetlib_except
        		messagefacility::MF_MessageLogger
			fhiclcpp::fhiclcpp
			${ROOT_BASIC_LIB_LIST}
        		ROOT::Geom
        		cetlib::cetlib
)
simple_plugin(
	      CRTGeometryHelper service
	      sbndcode_CRT
	      art::Framework_Core
	      art::Framework_IO_Sources
	      art::Framework_Principal
	      canvas::canvas
	      art::Framework_Services_Registry
	      messagefacility::MF_MessageLogger
)
