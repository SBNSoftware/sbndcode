set(
    MODULE_LIBRARIES
        sbncode::OpT0Finder_flashmatch_Base
        sbncode::OpT0Finder_flashmatch_Algorithms
        # sbndcode_OpDetReco_OpFlash_FlashFinder
        sbndcode_OpDetSim
        sbndcode::Geometry
        sbnobj::Common_Reco
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
        lardata::Utilities
        lardataobj::AnalysisBase
        lardataobj::RecoBase
        larsim::Simulation
        lardataobj::Simulation
        larpandora::LArPandoraInterface
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
        ROOT::Tree
)

cet_build_plugin(SBNDOpT0FinderAna art::module SOURCE SBNDOpT0FinderAna_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(SBNDOpT0Finder art::module SOURCE SBNDOpT0Finder_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()

add_subdirectory(job)
