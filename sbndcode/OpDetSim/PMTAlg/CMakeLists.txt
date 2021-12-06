art_make(
   TOOL_LIBRARIES
         art::Utilities
         messagefacility::MF_MessageLogger
         fhiclcpp::fhiclcpp
         cetlib::cetlib
         CLHEP::CLHEP
         cetlib::cetlib
         cetlib_except::cetlib_except
         
)


install_headers()
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
