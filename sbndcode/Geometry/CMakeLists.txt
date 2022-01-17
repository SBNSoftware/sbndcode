art_make(
          LIB_LIBRARIES larcorealg_Geometry
                        messagefacility::MF_MessageLogger
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        fhiclcpp::fhiclcpp
                        ${ROOT_BASIC_LIB_LIST}
                        ROOT::Geom
          SERVICE_LIBRARIES
                        sbndcode_Geometry
                        larcorealg_Geometry
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
