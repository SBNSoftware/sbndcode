set(
   TOOL_LIBRARIES
         art::Utilities
         larcore::Geometry_Geometry_service
         lardataobj::RawData
         lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
         sbndcode_Utilities_SignalShapingServiceSBND_service
         messagefacility::MF_MessageLogger
         fhiclcpp::fhiclcpp
         cetlib::cetlib
         CLHEP::CLHEP
         ROOT::Core
         art::Framework_Core
         art::Framework_Principal
         art::Framework_Services_Registry
         art::Persistency_Common
         art::Persistency_Provenance
         art::Utilities
         cetlib::cetlib
         cetlib_except::cetlib_except
)

cet_build_plugin(OpDeconvolutionAlgWiener art::tool SOURCE OpDeconvolutionAlgWiener_tool.cc LIBRARIES ${TOOL_LIBRARIES})

install_headers()
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
