set(
    MODULE_LIBRARIES
        sbndcode_OpDetReco_OpFlash_FlashFinder
        sbndcode_OpDetSim
        sbndcode_Geometry
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
        lardata::Utilities
        lardataobj::RecoBase
        larsim::Simulation
        lardataobj_Simulation
        nusimdata::SimulationBase
        canvas::canvas
        cetlib::cetlib
        cetlib_except::cetlib_except
        art::Framework_Core
        art::Framework_Principal
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
        ROOT::Core
) 

cet_build_plugin(SBNDFlashAna art::module SOURCE SBNDFlashAna_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(SBNDFlashFinder art::module SOURCE SBNDFlashFinder_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(SBNDMCFlash art::module SOURCE SBNDMCFlash_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()

add_subdirectory(FlashFinder)
add_subdirectory(FlashTools)
add_subdirectory(job)

