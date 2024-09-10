set(
          MODULE_LIBRARIES larcorealg::Geometry
                           larcore::Geometry_Geometry_service
                           larsim::Simulation
                           lardataobj::Simulation
                           larsim::MCCheater_BackTrackerService_service
                           larsim::MCCheater_ParticleInventoryService_service
                           lardata::Utilities
                           larevt::Filters
                           lardataobj::RawData
                           lardataobj::RecoBase
			   larreco::Calorimetry
                           larreco::RecoAlg
                           lardata::RecoObjects
                           larpandora::LArPandoraInterface
                           sbndcode_CRTUtils
                           sbndcode_CRT
                           sbnobj::Common_CRT
                           nusimdata::SimulationBase
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art_root_io::tfile_support
                           art_root_io::TFileService_service
                           art::Persistency_Common
                           art::Persistency_Provenance
                           messagefacility::MF_MessageLogger
                           fhiclcpp::fhiclcpp
                           ROOT::Geom
                           ROOT::XMLIO
                           ROOT::Gdml
                           ROOT::Core
                           sbndcode_RecoUtils
                           sbndcode_OpDetSim
			   sbndcode_CRT_CRTBackTracker
        )

cet_build_plugin(ToFAnalyzer art::module SOURCE ToFAnalyzer_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ToFFilter art::module SOURCE ToFFilter_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ToFProducer art::module SOURCE ToFProducer_module.cc LIBRARIES ${MODULE_LIBRARIES})

# install_headers()
install_fhicl()
install_source()

