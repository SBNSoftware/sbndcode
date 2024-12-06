add_subdirectory(tools)

set (MODULE_LIBRARIES
                        larcore::Geometry_Geometry_service
                        lardata::Utilities
                        lardata::ArtDataHelper
                        lardataobj::RawData
                        lardataobj::RecoBase
                        larevt::CalibrationDBI_Providers
                        sbndcode_Utilities_SignalShapingServiceSBND_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        ROOT::Core
                        ROOT::Tree
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
)
set (TOOL_LIBRARIES
                        larcore::Geometry_Geometry_service
                        lardata::Utilities
                        lardata::ArtDataHelper
                        lardataobj::RawData
                        lardataobj::RecoBase
                        sbndcode::Utilities_SignalShapingServiceSBND_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        ROOT::Core
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
)





cet_build_plugin(CalWireSBND art::module SOURCE CalWireSBND_module.cc LIBRARIES ${MODULE_LIBRARIES})

cet_build_plugin(HitEfficiencyAna art::module SOURCE HitEfficiencyAna_module.cc LIBRARIES ${MODULE_LIBRARIES})

cet_build_plugin(ROIFinderStandardSBND art::tool SOURCE ROIFinderStandardSBND_tool.cc LIBRARIES ${TOOL_LIBRARIES})

install_headers()
install_fhicl()
install_source()
