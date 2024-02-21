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

cet_build_plugin ( SBNDGeoHelper art::module
                     larcorealg::Geometry
                     larcore::Geometry_Geometry_service
                     art::Framework_Core
                     art::Framework_Principal
                     art::Framework_Services_Registry
                     art_root_io::tfile_support
                     ROOT::Core
                     art_root_io::TFileService_service
                     art::Framework_Services_Optional_RandomNumberGenerator_service
                     art::Persistency_Common
                     art::Persistency_Provenance
                     art::Utilities
		     canvas::canvas
                     messagefacility::MF_MessageLogger
                     fhiclcpp::fhiclcpp
                     cetlib::cetlib
	             cetlib_except::cetlib_except
                     ROOT::Geom
                     sbndcode_ChannelMaps_TPC_TPCChannelMapService_service

               BASENAME_ONLY
          )

install_headers()
install_fhicl()
install_source()
