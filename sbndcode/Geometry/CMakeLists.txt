art_make_library(
          SOURCE
                      ChannelMapSBNDAlg.cxx
                      GeoObjectSorterSBND.cxx
          LIBRARIES     larcorealg::Geometry
                        messagefacility::MF_MessageLogger
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        fhiclcpp::fhiclcpp
                        ROOT::Core
                        ROOT::Geom
)
cet_build_plugin( SBNDGeometryHelper art::service
          SOURCE SBNDGeometryHelper_service.cc
          LIBRARIES
                        sbndcode::Geometry
                        larcorealg::Geometry
                        larcore::Geometry_Geometry_service
                        art::Framework_Services_Registry
                        canvas::canvas
                        fhiclcpp::fhiclcpp
                        cetlib_except::cetlib_except
                        ROOT::Core
        )


install_headers()
install_fhicl()
install_source()

# GDML detector descriptions
add_subdirectory(gdml)
# Wrappers for easier access to geometry properties
add_subdirectory(GeometryWrappers)
