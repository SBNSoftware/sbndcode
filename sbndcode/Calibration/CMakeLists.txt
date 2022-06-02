add_subdirectory(tools)

art_make(MODULE_LIBRARIES
			larcore_Geometry_Geometry_service
                        lardata_Utilities
			lardata_ArtDataHelper
                        lardataobj_RawData
                        lardataobj_RecoBase
                        sbndcode_Utilities_SignalShapingServiceSBND_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
			art::Utilities			
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
                        ${ROOT_BASIC_LIB_LIST}

TOOL_LIBRARIES
			larcore_Geometry_Geometry_service
			lardata_Utilities
                        lardata_ArtDataHelper
                        lardataobj_RawData
                        lardataobj_RecoBase
                        sbndcode_Utilities_SignalShapingServiceSBND_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
			art::Utilities			
                        canvas::canvas
                        messagefacility::MF_MessageLogger

                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
                        ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()
