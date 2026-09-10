cet_enable_asserts()

set( nugraph_tool_lib_list
       lardataobj::RecoBase
       larrecodnn::NuGraphBaseTools
       TorchScatter::TorchScatter
       art::Framework_Core
       hep_concurrency::hep_concurrency
       art_plugin_types::tool
)

include_directories($ENV{HEP_HPC_INC})

simple_plugin(SBNDTrueNuSliceHitsProducer "module"
    LIBRARIES PRIVATE
    art::Framework_Core
    nusimdata::SimulationBase
    lardataobj::RecoBase
    lardataalg::DetectorInfo
    sbnobj::Common_Reco
    larsim::MCCheater_BackTrackerService_service           
    larsim::MCCheater_ParticleInventoryService_service
    lardata::DetectorClocksService
)

simple_plugin(SBNDHDF5Maker "module"
    LIBRARIES PRIVATE
    art::Framework_Core
    nusimdata::SimulationBase
    larsim::MCCheater_BackTrackerService_service           
    larsim::MCCheater_ParticleInventoryService_service
    hep_hpc_hdf5
    lardataobj::RecoBase
    lardataalg::DetectorInfo
    lardata::DetectorClocksService
    larcorealg::Geometry
    larcore::WireReadout
    ${HDF5_LIBRARIES}
)


install_headers()
install_fhicl()
install_source()
