################################################################################
###  CMakeLists.txt for gallery analysis executables
################################################################################
###
###  Usage:
###  1) edit the EXECUTABLE section adding executable names ("targets"),
###     the header paths they need, and the libraries they link to
###  2) if installing anywhere, choose where setting `GALLERY_INSTALL`
###     environment variable
###  3) the minimal setup for this example must include `larsoftobj` and
###     `sbndcode` (the latter may become optional depending on your code)
###  4) create a build directory and run `cmake` from there:
###      
###      mkdir -p build
###      cd build
###      cmake ..
###      make
###      
###
################################################################################
cmake_minimum_required(VERSION 3.10)

#
# making `cmake` aware of Clang may require some care;
# we pick `gallery` qualifiers and if they contain `cN` then it's Clang
# 
if ( $ENV{SETUP_GALLERY} MATCHES "-q +([^ ]+)( |$)" )
  set ( GalleryQualifiers ${CMAKE_MATCH_1} )
  if( ${GalleryQualifiers} MATCHES "(^|:)c[0-9]+(:|$)" )
    # make cmake acknowledge the difference between Apple Clang and standard Clang:
    cmake_policy(SET "CMP0025" "NEW")
    set(CMAKE_C_COMPILER "clang")
    set(CMAKE_CXX_COMPILER "clang++")
  else( ${GalleryQualifiers} MATCHES "(^|:)c[0-9]+(:|$)" )
    # don't bother, standard configuration will hopefully do
  endif( ${GalleryQualifiers} MATCHES "(^|:)c[0-9]+(:|$)" )
endif()

# make sure we know where to install things
if(DEFINED $ENV{GALLERY_INSTALL})
  set(GalleryInstallDir "$ENV{GALLERY_INSTALL}/" )
else(DEFINED $ENV{GALLERY_INSTALL})
  set(GalleryInstallDir "./" )
endif(DEFINED $ENV{GALLERY_INSTALL})

#
# declare the name of our project
#
project(galleryAnalysis CXX)


# install in the same directory as the source
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})

#
# include directories
#
include_directories($ENV{BOOST_INC})
include_directories($ENV{ROOT_INC} $ENV{CLHEP_INC})
include_directories($ENV{CETLIB_INC} $ENV{CETLIB_EXCEPT_INC} $ENV{CANVAS_INC} $ENV{CANVAS_ROOT_IO_INC})
include_directories($ENV{FHICLCPP_INC} $ENV{MESSAGEFACILITY_INC})
include_directories($ENV{HEP_CONCURRENCY_INC})
include_directories($ENV{GALLERY_INC})
include_directories($ENV{NUSIMDATA_INC} $ENV{LARCOREOBJ_INC} $ENV{LARDATAOBJ_INC})
include_directories($ENV{LARCOREALG_INC} $ENV{LARDATAALG_INC})

include_directories($ENV{SBNDCODE_INC})

#
# link directories and libraries
#
set(ROOT_LIBDIRS $ENV{ROOTSYS}/lib )
set(GALLERY_LIBDIRS ${ROOT_LIBDIRS} $ENV{MESSAGEFACILITY_LIB} $ENV{FHICLCPP_LIB} $ENV{CANVAS_LIB} $ENV{CANVAS_ROOT_IO_LIB} $ENV{CETLIB_EXCEPT_LIB} $ENV{CETLIB_LIB} $ENV{GALLERY_LIB} $ENV{HEP_CONCURRENCY_LIB})
set(LARSOFTOBJ_LIBDIRS $ENV{NUSIMDATA_LIB} $ENV{LARCOREOBJ_LIB} $ENV{LARDATAOBJ_LIB})
set(LARSOFT_LIBDIRS $ENV{CLHEP_LIB_DIR} $ENV{LARCOREALG_LIB} $ENV{LARDATAALG_LIB})
set(SBNDCODE_LIBDIRS $ENV{SBNDCODE_LIB})

link_directories(${GALLERY_LIBDIRS} ${LARSOFTOBJ_LIBDIRS} ${LARSOFT_LIBDIRS} ${SBNDCODE_LIBDIRS})

set(GALLERY_LIBS MF_MessageLogger fhiclcpp canvas cetlib_except cetlib gallery)
set(ROOT_CORELIBS Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread MultiProc pthread)
set(LARSOFTOBJ_LIBS nusimdata_SimulationBase larcoreobj_SummaryData lardataobj_RawData lardataobj_OpticalDetectorData lardataobj_RecoBase lardataobj_AnalysisBase lardataobj_MCBase lardataobj_Simulation)

set(SBNDGEOMETRY_LIBS larcorealg_Geometry sbndcode_Geometry)
set(SBNDDETINFO_LIBS lardataalg_DetectorInfo)

################################################################################
###  EXECUTABLE DIRECTIVES
################################################################################
add_executable(galleryAnalysis galleryAnalysis.cpp TrackAnalysis.cpp MCAssociations.cpp HitAnalysisAlg.cpp)
target_compile_features(galleryAnalysis PUBLIC cxx_std_17)
target_link_libraries(galleryAnalysis
  ${ROOT_CORELIBS} ${GALLERY_LIBS} ${LARSOFTOBJ_LIBS}
  ${SBNDGEOMETRY_LIBS} # if using geometry
  ${SBNDDETINFO_LIBS}  # if using LAr properties, detector clocks or detector properties
  )
install(TARGETS galleryAnalysis DESTINATION "${GalleryInstallDir}" )
