cet_enable_asserts()

set( sbnd_util_lib_list   lardata::Utilities_LArFFT_service
                          larcorealg::Geometry
                          larcore::Geometry_Geometry_service
                          lardata::Utilities
                          lardataalg::DetectorInfo
                          art::Framework_Core
                          art::Framework_Principal
                          art::Framework_Services_Registry
                          art::Persistency_Common
                          art::Persistency_Provenance 
                          art_root_io::tfile_support
                          art_root_io::TFileService_service
                          art::Framework_Services_System_TriggerNamesService_service
                          art::Utilities
                          canvas::canvas
                          messagefacility::MF_MessageLogger
                          cetlib::cetlib
                          cetlib_except::cetlib_except
                          ROOT::Geom
                          ROOT::Core
    )


cet_build_plugin( SignalShapingServiceSBND  art::service SOURCE SignalShapingServiceSBND_service.cc LIBRARIES
               ${sbnd_util_lib_list}
        ) 

install_headers()
install_fhicl()
install_source()
