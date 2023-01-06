set (MODULE_LIBRARIES
                larsim::Utils
                larsim::PhotonPropagation
                larsim::PhotonPropagation_PhotonVisibilityService_service
                larsim::LegacyLArG4
                larcorealg::GeoAlgo
                sbncode_OpT0Finder_flashmatch_Base
                larcorealg::Geometry
                larcore::Geometry_Geometry_service
                larsim::Simulation 
                lardataobj::Simulation
                larsim::MCCheater_BackTrackerService_service
                larsim::MCCheater_ParticleInventoryService_service
                lardata::Utilities
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
                art_root_io::tfile_support ROOT::Core
                art_root_io::TFileService_service
                art::Utilities canvas::canvas
                messagefacility::MF_MessageLogger
                
                fhiclcpp::fhiclcpp
                ROOT::Geom
                ROOT::XMLIO
                ROOT::Gdml
                ${ROOT_BASIC_LIB_LIST}
                sbndcode_RecoUtils
                sbndcode_OpDetSim
)
cet_build_plugin(LightCaloAna art::module SOURCE LightCaloAna_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()
