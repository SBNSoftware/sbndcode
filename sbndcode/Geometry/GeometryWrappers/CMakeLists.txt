art_make_library(LIBRARY_NAME TPCGeoAlg
  SOURCE TPCGeoAlg.cc
  LIBRARIES
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art::Framework_Principal
  messagefacility::MF_MessageLogger
  larcoreobj::geo_vectors
)

cet_build_plugin(CRTGeoService art::service
  SOURCE CRTGeoService_service.cc
  LIBRARIES
  messagefacility::MF_MessageLogger
  larcore::Geometry_Geometry_service
  larsim::Simulation
  nusimdata::SimulationBase
  sbndcode::CRTUtils
  sbndcode::ChannelMaps_CRT_CRTChannelMapService_service
  sbndcode::Calibration_CRT_CalibService_CRTCalibService_service
)

install_headers()
install_fhicl()
install_source()
