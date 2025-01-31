set(  
 MODULE_LIBRARIES    
 ${ART_FRAMEWORK_CORE}    
 ${ART_FRAMEWORK_PRINCIPAL}    
 ${ART_ROOT_IO_TFILE_SUPPORT} 
 ${ART_ROOT_IO_TFILESERVICE_SERVICE}     
 ${ART_FRAMEWORK_SERVICES_REGISTRY}    
 ${ART_ROOT_IO_TFILESERVICE_SERVICE}    
 ${FHICLCPP}    
 ${ROOT_BASIC_LIB_LIST}    
 ROOT::Tree
 ROOT::Core
 art_root_io::tfile_support
 art_root_io::TFileService_service
 art::Persistency_Common    
 art::Persistency_Provenance     
 art::Utilities    
 lardataobj::RecoBase
 larsim::MCCheater_BackTrackerService_service
 larcorealg::Geometry
 larcore::Geometry_Geometry_service
 larevt::SpaceCharge
 lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
 larsim::MCCheater_ParticleInventoryService_service
 sbnobj::SBND_CRT
 sbnobj::SBND_Timing
 sbndcode_GeoWrappers
 sbndcode_OpDetSim

 )

cet_build_plugin(FlashMatchAnalyzer art::Module SOURCE FlashMatchAnalyzer_module.cc LIBRARIES ${MODULE_LIBRARIES} )

install_source()
install_fhicl()
FILE(GLOB fcl_files *.fcl)

add_subdirectory(job)
add_subdirectory(EventIDFilter)
