
art_make(MODULE_LIBRARIES
    larreco_RecoAlg
    sbndcode_RecoUtils
    lardataobj_RecoBase
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ROOT::Core
    ROOT::Physics
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${MF_MESSAGELOGGER}
    ${ART_ROOT_IO_TFILE_SUPPORT}
)

install_headers()
install_fhicl()
install_source()
