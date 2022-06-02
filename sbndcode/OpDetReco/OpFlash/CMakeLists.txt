include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
art_make(
    LIB_LIBRARIES
        sbndcode_Geometry
        sbndcode_OpDetSim
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        lardata_Utilities
        lardataobj_RecoBase
        larsim_Simulation
        lardataobj_Simulation
        nusimdata::SimulationBase
        ${LARDATA_LIB}
        canvas::canvas
        cetlib::cetlib
        cetlib_except::cetlib_except
        ${ART_FRAMEWORK_ART}
        art::Framework_Core
        art::Framework_Principal
        ${ART_FRAMEWORK_UTILITIES}
        ${ART_FRAMEWORK_BASIC}
        art::Framework_Services_Registry
        art_root_io::tfile_support ROOT::Core
        art_root_io::TFileService_service
        art::Persistency_Common
        art::Persistency_Provenance
        art::Utilities
        messagefacility::MF_MessageLogger
        
        fhiclcpp::fhiclcpp
        ROOT::Geom
        ROOT::XMLIO
        ROOT::Gdml
        ${ROOT_BASIC_LIB_LIST}
    MODULE_LIBRARIES
        sbndcode_OpDetReco_OpFlash_FlashFinder
        sbndcode_OpDetSim
        sbndcode_Geometry
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        lardata_Utilities
        lardataobj_RecoBase
        larsim_Simulation
        lardataobj_Simulation
        nusimdata::SimulationBase
        ${LARDATA_LIB}
        canvas::canvas
        cetlib::cetlib
        cetlib_except::cetlib_except
        ${ART_FRAMEWORK_ART}
        art::Framework_Core
        art::Framework_Principal
        ${ART_FRAMEWORK_BASIC}
        art::Framework_Services_Registry
        art_root_io::tfile_support ROOT::Core
        art_root_io::TFileService_service
        art::Persistency_Common
        art::Persistency_Provenance
        art::Utilities
        messagefacility::MF_MessageLogger
        
        fhiclcpp::fhiclcpp
        ROOT::Geom
        ROOT::XMLIO
        ROOT::Gdml
        ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(FlashFinder)
add_subdirectory(job)

