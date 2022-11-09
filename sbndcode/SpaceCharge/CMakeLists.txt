art_make_library(
          SOURCE SpaceChargeSBND.cxx
          LIBRARIES
                        lardata::Utilities
			larcorealg::Geometry
                        messagefacility::MF_MessageLogger
                        larevt::SpaceCharge
                        cetlib::cetlib
                        cetlib_except::cetlib_except
			fhiclcpp::fhiclcpp
                        ROOT::Core
			
        ) 
install_headers()
install_fhicl()
install_source()
