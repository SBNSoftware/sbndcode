art_make(
          MODULE_LIBRARIES
                           larpandora::LArPandoraInterface
                           lardata::RecoObjects
                           lardataobj::RecoBase
                           lardataobj::AnalysisBase
                           larreco::RecoAlg
                           larcorealg::Geometry
                           larcore::Geometry_Geometry_service

                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art_root_io::tfile_support
                           art_root_io::TFileService_service
                           art::Persistency_Common
                           art::Persistency_Provenance

                           art::Utilities 
                           canvas::canvas

                           cetlib::cetlib
                           cetlib_except::cetlib_except
                           ${LARRECO_LIB}
                           ${LARDATA_LIB}
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_BASIC}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CLHEP}
                           sbndcode::RecoUtils
                           sbndcode::OpDetSim
                           sbnobj::Common_Reco
                            
                           messagefacility::MF_MessageLogger
                           fhiclcpp::fhiclcpp
                           ROOT::Geom
                           ROOT::XMLIO
                           ROOT::Gdml
                           ROOT::Core
                           ROOT::Tree

                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
                           ${ROOT_SPECTRUM}
                           ${ROOT_ROOFIT}
                           ${ROOT_ROOFITCORE}
        )
install_headers()
install_fhicl()
install_source()
