include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )
include_directories( $ENV{IFDHX_FQ_DIR}/inc )

add_subdirectory(CRTTools)
add_subdirectory(CRTUtils)
add_subdirectory(CRTAna)
add_subdirectory(CRTSimulation)
add_subdirectory(Examples)

art_make( 
    EXCLUDE
        CRTChannelMapAlg.cxx
        CRTGeometryHelper_service.cc
        CRTDetSim_module.cc
        CRTSimHitProducer_module.cc
        CRTTzeroProducer_module.cc
        CRTTrackProducer_module.cc
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

art_make_library( LIBRARY_NAME sbndcode_CRT
    SOURCE
        CRTChannelMapAlg.cxx
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

simple_plugin( CRTGeometryHelper service
        larcorealg_Geometry
        sbndcode_CRT
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        canvas::canvas
        cetlib_except::cetlib_except
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        
        ${ROOT_BASIC_LIB_LIST}
        cetlib::cetlib
)

simple_plugin( CRTDetSim module
        larcorealg_Geometry
        sbndcode_CRT
        sbndcode_CRT_CRTSimulation
        sbndcode_CRTData
        sbnobj_SBND_CRT
        sbnobj_Common_CRT
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        art::Persistency_Provenance
        canvas::canvas
        cetlib_except::cetlib_except
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
        messagefacility::MF_MessageLogger
        
        ${ROOT_BASIC_LIB_LIST}
        CLHEP::CLHEP
        cetlib::cetlib
)

simple_plugin( CRTSimHitProducer module
        larcorealg_Geometry
        sbndcode_CRT
        sbndcode_CRTData
        sbnobj_Common_CRT
        sbndcode_CRTUtils
        sbndcode_GeoWrappers
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
	art::Persistency_Provenance
        canvas::canvas
        cetlib_except::cetlib_except
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
	art_root_io::TFileService_service
        messagefacility::MF_MessageLogger
        ${ROOT_BASIC_LIB_LIST}
        CLHEP::CLHEP
        cetlib::cetlib
)

simple_plugin(CRTTrackProducer module
        larcorealg_Geometry
        sbndcode_CRT
        sbndcode_CRTData
        sbnobj_Common_CRT
        sbndcode_CRTUtils
        sbndcode_GeoWrappers
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
	art::Persistency_Provenance
        canvas::canvas
        cetlib_except::cetlib_except
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
	art_root_io::TFileService_service
        messagefacility::MF_MessageLogger
        
        ${ROOT_BASIC_LIB_LIST}
        CLHEP::CLHEP
        cetlib::cetlib
)

simple_plugin(CRTTzeroProducer module
        larcorealg_Geometry
        sbndcode_CRT
        sbndcode_CRTData
        sbnobj_Common_CRT
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
	art::Persistency_Provenance
        canvas::canvas
        cetlib_except::cetlib_except
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
	art_root_io::TFileService_service
        messagefacility::MF_MessageLogger
        
        ${ROOT_BASIC_LIB_LIST}
        CLHEP::CLHEP
        cetlib::cetlib
)

install_headers()
install_fhicl()
install_source()

