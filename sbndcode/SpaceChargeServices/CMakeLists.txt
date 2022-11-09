cet_build_plugin(SpaceChargeServiceSBND art::service SOURCE SpaceChargeServiceSBND_service.cc LIBRARIES
              sbndcode_SpaceCharge
	      art::Framework_Core
	      art::Framework_Principal
	      art::Framework_Services_Registry
	      art::Persistency_Common
	      art::Persistency_Provenance
	      art_root_io::tfile_support
	      art_root_io::TFileService_service
	      art::Utilities
) 
install_headers()
install_fhicl()
install_source()
