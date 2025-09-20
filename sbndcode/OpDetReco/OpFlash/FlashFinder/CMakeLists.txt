art_make_library(
    LIBRARY_NAME   sbndcode_OpDetReco_OpFlash_FlashFinder
    SOURCE
	FlashAlgoBase.cxx
	FlashAlgoFactory.cxx
	FlashFinderFMWKInterface.cxx
	FlashFinderManager.cxx
	PECalib.cxx
	SimpleFlashAlgo.cxx
     LIBRARIES
        sbndcode_Geometry
        sbndcode_OpDetSim
        larcore::Geometry_Geometry_service
        lardataobj::RecoBase
        larsim::Simulation
        lardataobj::Simulation
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

install_headers()
install_fhicl()
install_source()

