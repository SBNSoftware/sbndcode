art_make(
         MODULE_LIBRARIES
         ROOT::Tree
         ROOT::Graf3d
         CLHEP::Random
         art_root_io::TFileService_service
         art::Persistency_Common
         art::Persistency_Provenance
         nusimdata::SimulationBase
         nurandom::RandomUtils_NuRandomService_service
         larsim::MCCheater_ParticleInventoryService_service
         larsim::MCCheater_BackTrackerService_service
         larsim::Utils
         lardataobj::RecoBase
         lardata::Utilities
         sbnobj::Common_Reco
         sbndcode_PiZero_SecondShower
         sbndcode::GeoWrappers
)

install_fhicl()

add_subdirectory(SecondShower)
add_subdirectory(ShowerEnergyCorrection)
