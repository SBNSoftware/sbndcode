cet_enable_asserts()

set(
   MODULE_LIBRARIES
         sbndcode_OpDetSim
         larana::OpticalDetector_OpHitFinder
         larcore::Geometry_Geometry_service
         lardataobj::Simulation
         lardata::Utilities
         lardataobj::RawData
         lardataobj::RecoBase
         lardataobj::AnalysisBase
         lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
         larreco::Calorimetry
         sbndcode_Utilities_SignalShapingServiceSBND_service
         sbndcode_Utilities_SignalShapingServiceSBND_service
         nurandom::RandomUtils_NuRandomService_service
         art::Framework_Principal
         art::Framework_Services_Optional_RandomNumberGenerator_service
         canvas::canvas
         messagefacility::MF_MessageLogger
         fhiclcpp::fhiclcpp
         cetlib::cetlib
         cetlib_except::cetlib_except
         CLHEP::CLHEP
         ROOT::Core
) 

cet_build_plugin(SBNDOpHitFinder art::module SOURCE SBNDOpHitFinder_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()

add_subdirectory(job)
