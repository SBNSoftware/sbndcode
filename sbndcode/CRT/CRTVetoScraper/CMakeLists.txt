set(      MODULE_LIBRARIES larcorealg::Geometry
                           larcore::Geometry_Geometry_service
                           larsim::Simulation lardataobj::Simulation
                           larsim::MCCheater_BackTrackerService_service
                           larsim::MCCheater_ParticleInventoryService_service
 			   larsim::Utils
                           lardata::Utilities
                           larevt::Filters
                           lardataobj::RawData
                           lardataobj::RecoBase
			   larreco::Calorimetry
                           larreco::RecoAlg
                           lardata::RecoObjects
                           larpandora::LArPandoraInterface
                           nusimdata::SimulationBase
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art_root_io::tfile_support
                           art_root_io::TFileService_service
                           art::Persistency_Common
                           art::Persistency_Provenance
                           art::Utilities
                           messagefacility::MF_MessageLogger
                           ROOT::Core
                           ROOT::Tree
                           fhiclcpp::fhiclcpp
                           ROOT::Geom
                           ROOT::XMLIO
                           ROOT::Gdml
                           sbndcode_RecoUtils
                           sbndcode_OpDetSim

			   sbndaq_artdaq_core::sbndaq-artdaq-core_Obj_SBND
			   sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
			   sbnobj::SBND_Timing
    			   lardata::Utilities
    			   canvas::canvas
    			   lardataobj::headers
			   artdaq_core::artdaq-core_Data
                           sbnobj::SBND_CRT
                           sbnobj::SBND_Timing
                           sbndcode_GeoWrappers
                           sbndcode_CRT_CRTBackTracker


        )
cet_build_plugin(CRTVetoAnalysis art::Module SOURCE CRTVetoAnalysis_module.cc LIBRARIES ${MODULE_LIBRARIES} )
# install_headers()
install_fhicl()
install_source()


