set(	MODULE_LIBRARIES
 			sbnobj::Common_Trigger
            larcorealg::Geometry
			larcore::Geometry_Geometry_service
			lardataobj::RecoBase
			lardataobj::Simulation
			larcoreobj::SummaryData
			larana::OpticalDetector_OpHitFinder
			larreco::Calorimetry
			lardata::Utilities
			nusimdata::SimulationBase
			nurandom::RandomUtils_NuRandomService_service
			larsim::PhotonPropagation_PhotonVisibilityService_service
			art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art_root_io::tfile_support
			art_root_io::TFileService_service
			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
			canvas::canvas
			messagefacility::MF_MessageLogger
			messagefacility::headers
			fhiclcpp::fhiclcpp
			cetlib::cetlib
			CLHEP::Random
			ROOT::Geom
			ROOT::XMLIO
			ROOT::Gdml
			ROOT::FFTW
			ROOT::Core
            ROOT::Tree
)

cet_build_plugin(TPCPMTBarycenterMatchProducer art::module SOURCE TPCPMTBarycenterMatching_module.cc LIBRARIES ${MODULE_LIBRARIES})

add_subdirectory(job)

install_headers()
install_source()
install_fhicl()