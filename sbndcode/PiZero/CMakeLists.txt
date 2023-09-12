art_make(
         MODULE_LIBRARIES
         ROOT::Tree
         ROOT::Graf3d
         art_root_io::TFileService_service
         nusimdata::SimulationBase
         larsim::MCCheater_ParticleInventoryService_service
         larsim::MCCheater_BackTrackerService_service
         larsim::Utils
         lardataobj::RecoBase
         lardata::Utilities
         sbnobj::Common_Reco
)

install_fhicl()
