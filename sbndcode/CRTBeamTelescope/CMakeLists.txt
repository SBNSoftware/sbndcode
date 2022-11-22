add_subdirectory(fcls)
add_subdirectory(Geometry)
add_subdirectory(DataImport)
add_subdirectory(Reconstruction)

art_make(
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
                 sbnobj_Common_CRT
                 nusimdata::SimulationBase
                 art::Framework_Core
                 art::Framework_Principal
                 art::Framework_Services_Registry
                 art_root_io::tfile_support ROOT::Core
                 art_root_io::TFileService_service
                 art::Persistency_Common
                 art::Persistency_Provenance
                 art::Utilities
                 canvas::canvas
                 messagefacility::MF_MessageLogger

                 fhiclcpp::fhiclcpp
                 ROOT::Core
                 ROOT::Geom
                 ROOT::XMLIO
                 ROOT::Gdml
                 sbndcode_RecoUtils
                 sbndcode_OpDetSim
)
