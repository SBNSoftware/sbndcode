
set(TOOL_LIBRARIES
  larevt::CalibrationDBI_Providers
  lardataobj::RecoBase
  larreco::Calorimetry
  larcorealg::Geometry
  lardataalg::DetectorInfo
  lardata::ArtDataHelper
  canvas::canvas
  art::Framework_Services_Registry
  art::Utilities
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
  ROOT::Core
  ROOT::Hist
  wda::wda
)

cet_build_plugin(NormalizeDriftSQLite art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(NormalizeYZ art::tool LIBRARIES ${TOOL_LIBRARIES})

install_headers()
install_fhicl()
install_source()
