set(
      MODULE_LIBRARIES
                  larcorealg::Geometry
                  larcore::Geometry_Geometry_service
                  lardataobj::RecoBase
                  lardataobj::AnalysisBase 
                  lardataobj::RawData
                  lardata::RecoObjects 
                  lardata::Utilities
                  lardataobj::Simulation
                  nusimdata::SimulationBase
                  larcoreobj::SummaryData
                  art::Framework_Core
		  art::Framework_IO_Sources
		  art::Framework_Principal
		  art::Framework_Services_Registry
                  art_root_io::tfile_support ROOT::Core
		  art_root_io::TFileService_service
                  art::Persistency_Provenance
                  art::Persistency_Common
                  art::Utilities
                  canvas::canvas
		  ROOT::Geom
		  ROOT::Core
                  messagefacility::MF_MessageLogger
                  fhiclcpp::fhiclcpp
		  CLHEP::CLHEP
		  cetlib::cetlib
                  cetlib_except::cetlib_except
                  

)     

cet_build_plugin(HepMCFileGen art::module SOURCE HepMCFileGen_module.cc LIBRARIES ${MODULE_LIBRARIES})
install_headers()
install_source()
install_fhicl()
