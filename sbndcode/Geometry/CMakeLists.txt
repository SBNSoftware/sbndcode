art_make_library(
          SOURCE
                      GeoObjectSorterSBND.cxx
                      WireReadoutSorterSBND.cxx
          LIBRARIES     larcorealg::Geometry
                        messagefacility::MF_MessageLogger
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        fhiclcpp::fhiclcpp
                        ROOT::Core
                        ROOT::Geom
)

cet_build_plugin(GeoObjectSorterSBND art::tool LIBRARIES PRIVATE sbndcode::Geometry)
cet_build_plugin(WireReadoutSorterSBND art::tool LIBRARIES PRIVATE sbndcode::Geometry)

install_headers()
install_fhicl()
install_source()

# GDML detector descriptions
add_subdirectory(gdml)
# Wrappers for easier access to geometry properties
add_subdirectory(GeometryWrappers)
