
add_subdirectory(CRTTools)
add_subdirectory(CRTUtils)
#add_subdirectory(CRTAna)
add_subdirectory(CRTSimulation)
add_subdirectory(Examples)

art_make_library(
    SOURCE
        CRTGeoObjectSorter.cxx
    LIBRARY_NAME
        sbndcode_CRTData
    LIBRARIES
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
        lardataobj::RecoBase
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        canvas::canvas
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        fhiclcpp::fhiclcpp
        ROOT::Core
        ROOT::Geom
        cetlib::cetlib
)

art_make_library( LIBRARY_NAME sbndcode_CRT
    SOURCE
        CRTChannelMapAlg.cxx
    LIBRARIES
        larcorealg::Geometry
        sbndcode_CRTData
        cetlib_except::cetlib_except
        messagefacility::MF_MessageLogger
        fhiclcpp::fhiclcpp
        ROOT::Core
        ROOT::Geom
        cetlib::cetlib
)

cet_build_plugin( CRTGeometryHelper art::service SOURCE CRTGeometryHelper_service.cc
    LIBRARIES
        larcorealg::Geometry
        sbndcode_CRT
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        canvas::canvas
        cetlib_except::cetlib_except
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        ROOT::Core
        cetlib::cetlib
)

cet_build_plugin( CRTDetSim art::module SOURCE CRTDetSim_module.cc
    LIBRARIES
        larcorealg::Geometry
        sbndcode_CRT
        sbndcode_CRT_CRTSimulation
        sbndcode_CRTData
        sbnobj::SBND_CRT
        sbnobj::Common_CRT
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        art::Persistency_Provenance
        canvas::canvas
        cetlib_except::cetlib_except
        lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
        messagefacility::MF_MessageLogger
        ROOT::Core
        CLHEP::CLHEP
        cetlib::cetlib
)

cet_build_plugin( CRTSlimmer art::module SOURCE CRTSlimmer_module.cc
    LIBRARIES
        larcorealg::Geometry
        sbndcode_CRT
        sbndcode_CRT_CRTSimulation
        sbndcode_CRTData
        sbnobj::SBND_CRT
        sbnobj::Common_CRT
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        art::Persistency_Provenance
        canvas::canvas
        cetlib_except::cetlib_except
        lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
        messagefacility::MF_MessageLogger
        ROOT::Core
        CLHEP::CLHEP
        cetlib::cetlib
)

cet_build_plugin( CRTHitProducer art::module SOURCE CRTHitProducer_module.cc
    LIBRARIES
        larcorealg::Geometry
        sbndcode_CRT
        sbndcode_CRTData
        sbnobj::Common_CRT
        sbndcode_CRTUtils
        sbndcode_GeoWrappers
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
	art::Persistency_Provenance
        canvas::canvas
        cetlib_except::cetlib_except
        lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
	art_root_io::TFileService_service
        messagefacility::MF_MessageLogger
        ROOT::Core
        CLHEP::CLHEP
        cetlib::cetlib
)

cet_build_plugin(CRTTrackProducer art::module SOURCE CRTTrackProducer_module.cc
    LIBRARIES
        larcorealg::Geometry
        sbndcode_CRT
        sbndcode_CRTData
        sbnobj::Common_CRT
        sbndcode_CRTUtils
        sbndcode_GeoWrappers
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
	art::Persistency_Provenance
        canvas::canvas
        cetlib_except::cetlib_except
        lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
	art_root_io::TFileService_service
        messagefacility::MF_MessageLogger
        ROOT::Core
        CLHEP::CLHEP
        cetlib::cetlib
)

cet_build_plugin(CRTTzeroProducer art::module SOURCE CRTTzeroProducer_module.cc
    LIBRARIES
        larcorealg::Geometry
        sbndcode_CRT
        sbndcode_CRTData
        sbnobj::Common_CRT
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
	art::Persistency_Provenance
        canvas::canvas
        cetlib_except::cetlib_except
        lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
	art_root_io::TFileService_service
        messagefacility::MF_MessageLogger
        ROOT::Core
        CLHEP::CLHEP
        cetlib::cetlib
)

simple_plugin(CRTEventDisplay module
        larcorealg::Geometry
        sbndcode_CRT
        sbndcode_CRTData
        sbndcode_CRTUtils
        sbnobj_Common_CRT
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
	art::Persistency_Provenance
        canvas::canvas
        cetlib_except::cetlib_except
        lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
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

