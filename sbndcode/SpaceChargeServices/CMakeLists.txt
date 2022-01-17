simple_plugin(SpaceChargeServiceSBND "service"
              sbndcode_SpaceCharge
	      art::Framework_Core
	      art::Framework_Principal
	      art::Framework_Services_Registry
	      art::Persistency_Common canvas
	      art::Persistency_Provenance canvas
	      ${ART_FRAMEWORK_SERVICES_BASIC}
	      art_root_io::tfile_support ROOT::Core
	      art_root_io::TFileService_service
	      ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
	      art::Utilities canvas
	      messagefacility::MF_MessageLogger
	      
  	      cetlib::cetlib cetlib_except
   	      ${ROOT_BASIC_LIB_LIST}
)
install_headers()
install_fhicl()
install_source()
