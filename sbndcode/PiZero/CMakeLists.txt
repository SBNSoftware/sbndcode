art_make(
	 MODULE_LIBRARIES
	 ROOT::Tree
	 art_root_io::TFileService_service
	 nusimdata::SimulationBase
	 larsim::MCCheater_ParticleInventoryService_service
	 larsim::MCCheater_BackTrackerService_service
	 larsim::Utils
	 lardataobj::RecoBase
	 lardata::Utilities
)

install_fhicl()
