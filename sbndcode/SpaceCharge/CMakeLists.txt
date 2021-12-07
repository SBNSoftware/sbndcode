art_make( NO_PLUGINS          
          LIB_LIBRARIES
                        lardata_Utilities
			larcorealg_Geometry
                        messagefacility::MF_MessageLogger
                        
                        cetlib::cetlib cetlib_except
			fhiclcpp::fhiclcpp
                        ${ROOT_BASIC_LIB_LIST}
			
        )
install_headers()
install_fhicl()
install_source()
