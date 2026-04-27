# Data-driven cathode volume library
art_make_library(
  LIBRARY_NAME sbndcode_CathodeFilter
  SOURCE       SBNDCathode.cxx
  LIBRARIES    
   cetlib::cetlib
   lardataobj::Simulation
   art_root_io::TFileService_service
   ROOT::Core
   ROOT::GenVector
  INTERFACE
)

# The art producer module
cet_build_plugin(SBNDCathode  art::tool
  SOURCE 
    SBNDCathode.cxx
  LIBRARIES
    sbndcode_CathodeFilter
    lardataobj::Simulation
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art::Persistency_Common
    art::Utilities
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    CLHEP::CLHEP
    ROOT::Core
    ROOT::GenVector
)


# The art producer module
cet_build_plugin( CathodeFilterSBND art::module
  SOURCE
    CathodeSimEnergyDepositFilter_module.cc
  LIBRARIES
    sbndcode_CathodeFilter
    larcorealg::Geometry
    larcore::Geometry_Geometry_service
    larsim::Simulation
    lardataobj::Simulation
    larcoreobj::headers
    lardata::Utilities
    lardata::RecoObjects
    nusimdata::SimulationBase
    lardataobj::Simulation
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art::Persistency_Common
    art::Utilities
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    CLHEP::CLHEP
    ROOT::Core
    art_root_io::tfile_support
    art_root_io::TFileService_service
    ROOT::MathCore
    ROOT::GenVector    
)

#cet_build_plugin(CathodeFilterSBND art::module SOURCE CathodeSimEnergyDepositFilter_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_source()
install_fhicl(SUBDIRS fcl)
cet_enable_asserts()

# install sbnd_pds_mapping.json with mapping of the photon detectors                                                                                    
install_fw(LIST sbnd_cathode_v6.txt)