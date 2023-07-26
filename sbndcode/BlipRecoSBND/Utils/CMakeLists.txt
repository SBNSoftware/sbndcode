art_make_library(LIBRARY_NAME  sbndcode_BlipUtils
                 SOURCE    BlipUtils.cc
                 LIBRARIES larcorealg::Geometry
                           larcore::Geometry_Geometry_service
                           larsim::Simulation
                           lardataobj::Simulation
                           larsim::MCCheater_BackTrackerService_service
                           larsim::MCCheater_ParticleInventoryService_service
                           lardata::Utilities
                           larevt::Filters
                           lardataobj::RawData
			                     larreco::Calorimetry
                           lardataobj::RecoBase
                           lardata::RecoObjects
                           larpandora::LArPandoraInterface
                           nusimdata::SimulationBase
                           cetlib::cetlib
                           cetlib_except::cetlib_except
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art_root_io::tfile_support
                           art_root_io::TFileService_service
                           art::Persistency_Common
                           canvas::canvas
                           art::Persistency_Provenance
                           art::Utilities
                           messagefacility::MF_MessageLogger
                           ROOT::Core
                           fhiclcpp::fhiclcpp
                           ROOT::Geom
                           ROOT::XMLIO
                           ROOT::Gdml
                           sbndcode_CRT
                           sbndcode_CRTUtils
                           sbnobj::Common_CRT
                           sbndcode_CosmicIdUtils
)

install_headers()
install_source()
install_fhicl()
