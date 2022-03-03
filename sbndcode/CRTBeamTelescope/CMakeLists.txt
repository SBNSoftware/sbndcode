add_subdirectory(fcls)
add_subdirectory(Geometry)

art_make(
MODULE_LIBRARIES larcorealg_Geometry
                 larcore_Geometry_Geometry_service
                 larsim_Simulation lardataobj_Simulation
                 larsim_MCCheater_BackTrackerService_service
                 larsim_MCCheater_ParticleInventoryService_service
                 lardata_Utilities
                 larevt_Filters
                 lardataobj_RawData
                 lardataobj_RecoBase
                 larreco_Calorimetry
                 larreco_RecoAlg
                 lardata_RecoObjects
                 larpandora_LArPandoraInterface
                 sbndcode_CRTUtils
                 sbndcode_CRT
                 sbnobj_Common_CRT
                 nusimdata::SimulationBase
                 art::Framework_Core
                 art::Framework_Principal
                 art::Framework_Services_Registry
                 art_root_io::tfile_support ROOT::Core
                 art_root_io::TFileService_service
                 art::Persistency_Common canvas
                 art::Persistency_Provenance canvas
                 art::Utilities canvas
                 messagefacility::MF_MessageLogger

                 fhiclcpp::fhiclcpp
                 ROOT::Geom
                 ROOT::XMLIO
                 ROOT::Gdml
                 ${ROOT_BASIC_LIB_LIST}
                 sbndcode_RecoUtils
                 sbndcode_OpDetSim
)

