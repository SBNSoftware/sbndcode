cet_build_plugin(CherenkovPhotonConverter art::module
      LIBRARIES PRIVATE
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
        lardataobj::Simulation
        artg4tk::pluginDetectors_gdml
        sbndcode::OpDetSim
)

install_source()
install_fhicl()
