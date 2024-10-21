
art_make(
    MODULE_LIBRARIES sbnobj::Common_SBNEventWeight
                     larcorealg::CoreUtils
                     larcorealg::Geometry
                     larcore::Geometry_Geometry_service
                     larsim::Simulation
                     lardataobj::Simulation
                     lardata::Utilities
                     nusimdata::SimulationBase
                     art::Framework_Core
                     art::Framework_Principal
                     art::Framework_Services_Registry
                     art_root_io::tfile_support
                     art_root_io::TFileService_service
                     art::Persistency_Common
                     art::Persistency_Provenance
                     art::Utilities
                     messagefacility::MF_MessageLogger
                     fhiclcpp::fhiclcpp
                     ROOT::Core
                     ROOT::Tree
                     ROOT::Geom
                     ROOT::XMLIO
                     ROOT::Gdml
                     ROOT::EG
)

# install_headers()
install_fhicl()
install_source()
