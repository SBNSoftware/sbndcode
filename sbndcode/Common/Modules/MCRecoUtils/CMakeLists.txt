art_make(
             LIB_LIBRARIES larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation lardataobj_Simulation
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           lardata_Utilities
                           lardataobj_RawData
                           lardataobj_RecoBase
			   lardata_RecoObjects  
                           nusimdata_SimulationBase
                           cetlib cetlib_except
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           art_Persistency_Common canvas
                           art_Persistency_Provenance canvas
                           art_Utilities canvas
                           ${MF_MESSAGELOGGER}
                           ${FHICLCPP}
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
                           ${Boost_SYSTEM_LIBRARY}
	                   ${ART_ROOT_IO_TFILE_SUPPORT}	
			   ${ART_ROOT_IO_TFILESERVICE_SERVICE}


        )

install_headers()
install_fhicl()
install_source()


