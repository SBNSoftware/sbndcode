art_make(
   TOOL_LIBRARIES
         art_Utilities
         ${MF_MESSAGELOGGER}
         ${FHICLCPP}
         ${CETLIB}
         ${CLHEP}
         cetlib
         cetlib_except
         ${MF_UTILITIES}
)


install_headers()
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
