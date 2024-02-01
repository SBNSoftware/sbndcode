art_make(
         MODULE_LIBRARIES
         ROOT::Tree
         art_root_io::TFileService_service
         art::Persistency_Common
         art::Persistency_Provenance
         nusimdata::SimulationBase
         larsim::MCCheater_ParticleInventoryService_service
         larsim::Utils
         lardataobj::RecoBase
         lardata::Utilities
)

install_fhicl()
