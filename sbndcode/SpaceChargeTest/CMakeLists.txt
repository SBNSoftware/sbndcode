simple_plugin(SpaceChargeTest "module"
                           larcore_Geometry_Geometry_service
                           lardata_Utilities
                           larevt_Filters
                           lardataobj_RawData
                           lardataobj_RecoBase
                           lardata_RecoObjects  
                           larpandora_LArPandoraInterface
                           nusimdata::SimulationBase
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art_root_io::tfile_support ROOT::Core
                           art_root_io::TFileService_service
                           art::Persistency_Common canvas
                           art::Persistency_Provenance canvas
                           art::Utilities canvas
                           messagefacility::MF_MessageLogger
                           
                           fhiclcpp::fhiclcpp
                           ROOT::Geom
                           ROOT::XMLIO
                           ROOT::Gdml
                           ${ROOT_BASIC_LIB_LIST}
                           sbndcode_RecoUtils
			   sbndcode_SpaceCharge
        )

install_headers()
install_fhicl()
install_source()
