set (MODULE_LIBRARIES
                larsim::PhotonPropagation
                larsim::OpticalPath
                larsim::LegacyLArG4
                larcorealg::Geometry
                larcore::Geometry_Geometry_service
                lardata::Utilities
                lardataobj::RecoBase
                art::Framework_Core
                art::Framework_Principal
                art::Framework_Services_Registry
                art::Utilities canvas::canvas
                art_root_io::TFileService_service
                art_root_io::tfile_support 
                messagefacility::MF_MessageLogger
                sbnobj::Common_Reco
                fhiclcpp::fhiclcpp
                sbndcode_OpDetSim
                ROOT::Core 
                ROOT::Tree
)

cet_build_plugin(LightCaloProducer art::module SOURCE LightCaloProducer_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()
