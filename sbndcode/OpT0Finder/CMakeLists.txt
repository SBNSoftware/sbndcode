# include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# Turning off WERROR only for this package, because it needs to
# call OpFlastScintillation which has an unused variable.
# Once the semi-analytic photon propogation will be moved
# to its own service, and we won't need to call OpFlastScintillation
# in larsim, we will be able to remobe this. -- Marco Del Tutto, Oct 2020
# cet_set_compiler_flags(DIAGS CAUTIOUS
#   # WERROR
#   NO_UNDEFINED
#   EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
# )
include_directories($ENV{GEANT4_FQ_DIR}/include)

art_make(
    LIB_LIBRARIES
        sbncode_OpT0Finder_flashmatch_Base
        sbncode_OpT0Finder_flashmatch_Algorithms
        sbnbcode_Geometry
        sbnbcode_OpDetSim
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        lardata_Utilities
        lardataobj_RecoBase
        larsim_Simulation
        # larsim_LegacyLArG4
        lardataobj_AnalysisBase
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
        sbncode_OpT0Finder_flashmatch_Base
        sbncode_OpT0Finder_flashmatch_Algorithms
        # sbndcode_OpDetReco_OpFlash_FlashFinder
        sbndcode_OpDetSim
        sbndcode_Geometry
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        lardata_Utilities
        lardataobj_AnalysisBase
        lardataobj_RecoBase
        larsim_Simulation
        # larsim_LegacyLArG4
        lardataobj_Simulation
        larpandora_LArPandoraInterface
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

add_subdirectory(job)
