# make a shared library from OnlinePMTAna.cxx
cet_make_library( LIBRARY_NAME SBNOM_OnlinePMTAna
  SOURCE  
    OnlinePMTAna.cxx
  LIBRARIES
    ${LARDATAOBJ}
    lardataobj_RawData
    ${ART_UTILITIES}
    ${FHICLCPP}
    ${ART_FRAMEWORK_CORE}
    ${ROOT_CORE}
    ${ROOT_BASIC_LIB_LIST}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_BASIC}
    ${ART_ROOT_IO_TFILE_SUPPORT} 
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_PERSISTENCY_PROVENANCE}
    canvas
    cetlib cetlib_except
    art_Persistency_Common
)

# make a module from OnlinePMTAna_module.cc
simple_plugin( OnlinePMTAna module
  SBNOM_OnlinePMTAna
  ${LARDATAOBJ}
  lardataobj_RawData
  ${ART_UTILITIES}
  ${FHICLCPP}
  ${ART_FRAMEWORK_CORE}
  ${ROOT_BASIC_LIB_LIST}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  art_Persistency_Common
)

