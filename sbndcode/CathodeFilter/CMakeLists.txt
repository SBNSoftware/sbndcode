include_directories( $ENV{LARDATAOBJ_INC} )

# Data-driven cathode volume library
cet_make_library(
  LIBRARY_NAME sbndcode_CathodeFilter
  SOURCE       SBNDCathode.cxx
  LIBRARIES    INTERFACE
)

# The art producer module
cet_build_plugin( CathodeSimEnergyDepositFilter art::module
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
)

install_headers()
install_source()
install_fhicl(SUBDIRS fcl)

# install the data file under fw-search-path so FindFile can locate it
add_subdirectory(data)
