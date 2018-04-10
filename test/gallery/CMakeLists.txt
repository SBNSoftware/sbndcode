# for PREBUILT tests, the test name is the executable name, and no argument is allowed

# Only for GCC compiler (eXX qualifiers)
if($ENV{MRB_QUALS} MATCHES "^(.*:)?e[0-9]+(:.*)?$" )
  cet_test(compileGalleryAnalysis.sh PREBUILT)
endif()
