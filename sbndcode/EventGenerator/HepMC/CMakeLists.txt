include_directories( $ENV{IFDHC_FQ_DIR}/inc )

set( IFDH_AS_SERVICE_LIB_LIST
     art_Framework_Services_Registry
     ${Boost_SYSTEM_LIBRARY}
     ${IFDH_SERVICE} )

include_directories ( . )

art_make(
      LIBRARY_NAME     
                  sbndcode_hepmc
      LIB_LIBRARIES
                  larcorealg_Geometry
                  larcore_Geometry_Geometry_service
                  lardataobj_RecoBase
                  lardataobj_AnalysisBase 
                  lardataobj_RawData
                  lardata_RecoObjects 
                  lardata_Utilities
                  lardataobj_Simulation
                  nusimdata_SimulationBase
                  larcoreobj_SummaryData
                  lardataobj_RawData
                  lardataobj_RecoBase
                  lardataobj_AnalysisBase
                  ${ART_FRAMEWORK_CORE}
                  ${ART_FRAMEWORK_BASIC}
                  ${ART_FRAMEWORK_IO_PRODUCTMIX}
                  ${ART_FRAMEWORK_IO_ROOT}
                  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                  ${ART_FRAMEWORK_MODULES}
                  ${ART_FRAMEWORK_PRINCIPAL}
                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                  ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                  art_Persistency_Provenance
                  art_Persistency_Common
                  art_Utilities
                  canvas
                  ifdh
                  ${IFDH}
                  ${IFDH_AS_SERVICE_LIB_LIST}
                  ${MF_MESSAGELOGGER}
                  ${MF_UTILITIES}
                  ${BOOST_LIB}
                  ${Boost_SYSTEM_LIBRARY}
                  ${Boost_FILESYSTEM_LIBRARY}
                  ${FHICLCPP}
                  ${CETLIB}
                  cetlib cetlib_except
                  ${ROOT_BASIC_LIB_LIST}
                  

      MODULE_LIBRARIES
                  larcorealg_Geometry
                  larcore_Geometry_Geometry_service
                  lardataobj_RecoBase
                  lardataobj_AnalysisBase 
                  lardataobj_RawData
                  lardata_RecoObjects 
                  lardata_Utilities
                  lardataobj_Simulation
                  nusimdata_SimulationBase
                  larcoreobj_SummaryData
                  lardataobj_RawData    
                  lardataobj_RecoBase
                  lardataobj_AnalysisBase
                  ${ART_FRAMEWORK_CORE}
                  ${ART_FRAMEWORK_BASIC}
                  ${ART_FRAMEWORK_IO_PRODUCTMIX}
                  ${ART_FRAMEWORK_IO_ROOT}
                  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                  ${ART_FRAMEWORK_MODULES}
                  ${ART_FRAMEWORK_PRINCIPAL}
                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                  ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                  art_Persistency_Provenance
                  art_Persistency_Common
                  art_Utilities
                  canvas
                  ifdh
                  ${IFDH}
                  ${IFDH_AS_SERVICE_LIB_LIST}
                  ${MF_MESSAGELOGGER}
                  ${MF_UTILITIES}
                  ${BOOST_LIB}
                  ${Boost_SYSTEM_LIBRARY}
                  ${Boost_FILESYSTEM_LIBRARY}
                  ${FHICLCPP}
                  ${CETLIB}
                  cetlib cetlib_except
                  ${ROOT_BASIC_LIB_LIST}

)     

install_headers()
install_source()
install_fhicl()
