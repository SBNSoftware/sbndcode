art_make(
         MODULE_LIBRARIES
         ROOT::Tree
         ROOT::Graf3d
         art_root_io::TFileService_service
         art::Persistency_Common
         art::Persistency_Provenance
         nusimdata::SimulationBase
         larsim::MCCheater_ParticleInventoryService_service
         larsim::MCCheater_BackTrackerService_service
         larsim::Utils
         lardataobj::RecoBase
         lardata::Utilities
)

install_fhicl()
