include (phot::OpticalPath)

cet_build_plugin(SBNDOpticalPath phot::OpticalPath
  LIBRARIES PRIVATE
  larsim::OpticalPath
  larcoreobj::geo_vectors
)

install_headers()
install_fhicl()
install_source()