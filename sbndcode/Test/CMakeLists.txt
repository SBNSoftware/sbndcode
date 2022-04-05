art_make(MODULE_LIBRARIES
    larreco_RecoAlg
    lardataobj_RecoBase
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ROOT::Core
    ROOT::Physics
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${MF_MESSAGELOGGER}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    sbnobj_Common_CRT
)

install_headers()
install_fhicl()
install_source()
