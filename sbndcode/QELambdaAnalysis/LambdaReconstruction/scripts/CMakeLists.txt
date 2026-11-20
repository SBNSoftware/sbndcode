# using/not using larsoft
if(DEFINED ENV{LARSOFT_VERSION})
    message(status " LARSOFT DEFINED" )
    file(GLOB xml_file *.xml )
    install_fw(LIST ${xml_file})
endif()
