cet_enable_asserts()

set( sbnd_util_lib_list lardata_Utilities_LArFFT_service
                          larcorealg_Geometry
                          larcore_Geometry_Geometry_service
                          lardata_Utilities
                          lardataalg_DetectorInfo
                          art::Framework_Core
                          art::Framework_Principal
                          art::Framework_Services_Registry
                          art::Persistency_Common
                          art::Persistency_Provenance 
                          ${ART_FRAMEWORK_SERVICES_BASIC}
                          art_root_io::tfile_support ROOT::Core
                          art_root_io::TFileService_service
                          art_Framework_Services_System_TriggerNamesService_service
                          art::Utilities
                          canvas::canvas
                          messagefacility::MF_MessageLogger
                          cetlib::cetlib cetlib_except
                          ${ROOT_BASIC_LIB_LIST}
                          ROOT::Geom
    )


simple_plugin( SignalShapingServiceSBND  "service"
               ${sbnd_util_lib_list}
        )

install_headers()
install_fhicl()
install_source()
