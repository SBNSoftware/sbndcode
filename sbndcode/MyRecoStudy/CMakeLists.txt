set(MODULE_LIBRARIES
          sbnci_Common_Modules_MCRecoUtils
          lardataobj::RecoBase
          larsim::Utils
	  lardataalg::DetectorInfo
          nusimdata::SimulationBase
          ROOT::Core
          ROOT::Physics
          ROOT::Tree
          messagefacility::MF_MessageLogger
          art::Framework_Services_Registry
          art_root_io::tfile_support
          art_root_io::TFileService_service
        )

cet_build_plugin(PFPValidation art::Module SOURCE PFPValidation_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_fhicl()
install_source()
