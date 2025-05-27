set(MODULE_LIBRARIES
	larsim::Utils
	lardataobj::RecoBase
    lardataobj::MCBase
    nusimdata::SimulationBase
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_ROOT_IO_TFILE_SUPPORT} 
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::tfile_support
    art_root_io::TFileService_service
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities

    ${ROOT_BASIC_LIB_LIST}

	sbndcode_GeoWrappers
	sbndcode_RecoUtils

    larcorealg::Geometry
    larcore::Geometry_Geometry_service
    lardata::Utilities
    larevt::Filters
    lardataobj::RawData
    lardataobj::RecoBase
    lardataobj::MCBase
    larcoreobj::SummaryData
    larsim::MCCheater_BackTrackerService_service
    larsim::MCCheater_ParticleInventoryService_service
    nusimdata::SimulationBase
    art::Persistency_Common canvas::canvas
    art::Persistency_Provenance
    art::Utilities
    cetlib::cetlib cetlib_except::cetlib_except
    ROOT::Geom
    ROOT::Core
    ROOT::Tree
	)

cet_build_plugin(RecoAnalysis art::module SOURCE RecoAnalysis_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(RecoFilter art::module SOURCE RecoFilter_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(AlphaFilter art::module SOURCE AlphaFilter_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()