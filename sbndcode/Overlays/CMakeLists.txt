cet_enable_asserts()

cet_build_plugin(OverlayProducts art::module LIBRARIES
        #icaruscode::Timing
        sbnobj::Common_PMT_Data
        larcore::Geometry_Geometry_service
        art_root_io::TFileService_service
        art_root_io::tfile_support
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        fhiclcpp::fhiclcpp
        cetlib_except::cetlib_except
        ROOT::Tree
        lardata::Utilities
        lardataobj::RecoBase
        lardata::ArtDataHelper
)

cet_build_plugin(OverlayAna
        art::module LIBRARIES
        sbnobj::Common_PMT_Data
        larcore::Geometry_Geometry_service
        art_root_io::TFileService_service
        art_root_io::tfile_support
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        fhiclcpp::fhiclcpp
        cetlib_except::cetlib_except
        ROOT::Tree
        lardata::Utilities
        lardataobj::RecoBase
        lardata::ArtDataHelper
)

install_headers()
install_fhicl()
install_source()

