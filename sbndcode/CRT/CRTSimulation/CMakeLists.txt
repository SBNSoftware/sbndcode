art_make(   LIB_LIBRARIES  larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           lardataobj_Simulation
                           lardataobj_RawData
                           lardataobj_RecoBase
                           lardata_RecoObjects
                           lardataalg_DetectorInfo
                           cetlib::cetlib cetlib_except
                           messagefacility::MF_MessageLogger
                           fhiclcpp::fhiclcpp

                           ${ROOT_BASIC_LIB_LIST}
                           ${ROOT_MATHMORE}

                           sbnobj_SBND_CRT
                           sbndcode_GeoWrappers
        )

install_headers()
install_fhicl()
install_source()
