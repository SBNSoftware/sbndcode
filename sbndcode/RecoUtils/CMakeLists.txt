art_make(
             LIB_LIBRARIES larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation lardataobj_Simulation
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           lardata_Utilities
                           larevt_Filters
                           lardataobj_RawData
                           lardataobj_RecoBase
                           lardata_RecoObjects  
                           larpandora_LArPandoraInterface
                           nusimdata::SimulationBase
                           cetlib::cetlib cetlib_except
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art_root_io::tfile_support ROOT::Core
                           art_root_io::TFileService_service
                           art::Persistency_Common canvas
                           art::Persistency_Provenance canvas
                           art::Utilities canvas
                           messagefacility::MF_MessageLogger
                           
                           fhiclcpp::fhiclcpp
                           ROOT::Geom
                           ROOT::XMLIO
                           ROOT::Gdml
                           ${ROOT_BASIC_LIB_LIST}
                           
        )

install_headers()
install_fhicl()
install_source()


#art_make( LIBRARY_NAME RecoUtils
#	  LIB_LIBRARIES larcorealg_Geometry
#	  	        larcore_Geometry_Geometry_service
#		        larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
# 	  	        larsim_MCCheater_BackTrackerService_service
#	  	        lardata_Utilities
#larevt_Filters
#	  	        lardataobj_RawData
#	  	        lardataobj_RecoBase
#	  	        lardata_RecoObjects
#larreco_RecoAlg
#larreco_Calorimetry
#larpandora_LArPandoraInterface
#          	        nusimdata::SimulationBase
#          	        ${LARRECO_LIB}
#          	        ${LARDATA_LIB}
#	  	        messagefacility::MF_MessageLogger
#	  	        
#	  	        fhiclcpp::fhiclcpp
#	  	        cetlib::cetlib cetlib_except
#          	        ROOT::Geom
#          	        ROOT::XMLIO
#          	        ROOT::Gdml
#	  	        ${ROOT_BASIC_LIB_LIST}
#	  MODULE_LIBRARIES RecoUtils
#          	           art::Framework_Core
#	  	           art::Framework_Principal
#		       	   art::Framework_Services_Registry
#	  	       	   art_root_io::tfile_support ROOT::Core
#			   art_root_io::TFileService_service
#	  	           art::Persistency_Common canvas
#	  	           art::Persistency_Provenance canvas
#	  	           art::Utilities canvas
#        )
#
#install_headers()
##install_fhicl()
#install_source()
##install_scripts()
