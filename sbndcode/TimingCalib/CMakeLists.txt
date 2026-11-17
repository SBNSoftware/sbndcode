# dk2nu::Tree (MC timing offset); sbndcode's top-level CMakeLists does not find it
find_package( dk2nudata REQUIRED )

set (
  MODULE_LIBRARIES
    lardataobj::RecoBase
    lardataobj::Simulation
    lardata::DetectorClocksService
    larsim::MCCheater_BackTrackerService_service
    larsim::MCCheater_ParticleInventoryService_service
    nusimdata::SimulationBase
    larcoreobj::SummaryData
    dk2nu::Tree
    sbnobj::Common_POTAccounting
    sbnobj::Common_Reco
    sbnobj::SBND_Timing
    art_root_io::TFileService_service
    messagefacility::MF_MessageLogger
    cetlib_except::cetlib_except
    ROOT::Tree
)

cet_build_plugin(TimingCalibAna art::EDAnalyzer LIBRARIES PRIVATE ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()
