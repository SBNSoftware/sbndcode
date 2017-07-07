# include_directories ( $ENV{XERCES_C_INC} )

art_make_exec(NAME preparse
  SOURCE preparse.cpp
  LIBRARIES
    ${ROOT_BASIC_LIB_LIST}
  )

install_source()
