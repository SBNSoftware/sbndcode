add_subdirectory(CRTAna)
add_subdirectory(CRTBackTracker)
add_subdirectory(CRTEventDisplay)
add_subdirectory(CRTReco)
add_subdirectory(CRTSimulation)
add_subdirectory(CRTTPCMatching)
add_subdirectory(CRTUtils)

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

cet_build_plugin( CRTAuxDetInitializer art::tool
    LIBRARIES REG
        larcorealg::Geometry
        sbndcode_CRTData
        cetlib_except::cetlib_except
        messagefacility::MF_MessageLogger
        fhiclcpp::fhiclcpp
        cetlib::cetlib
)

install_headers()
install_fhicl()
install_source()
