################################################################################
###  CMakeLists.txt for gallery analysis executables
################################################################################
###
###  Usage:
###  1) edit the EXECUTABLE section adding executable names ("targets"),
###     the header paths they need, and the libraries they link to
###  2) the minimal setup for this example must include `larsoftobj` and
###     `sbndcode` (the latter may become optional depending on your code)
###  3) create a build directory and run `cmake` from there:
###      
###      mkdir -p build
###      cd build
###      cmake ..
###      make
###      
###     To build in debug mode add `-DCMAKE_BUILD_TYPE=Debug` to `cmake` command
###     above (`Release` sets the optimised version).
###
###  The required development setup includes:
###   * UPS databases (with canvas, sbndcode etc.)
###   * cetmodules (required to find some UPS libraries; pulls in cmake)
###   * the UPS products mentioned above (`sbndcode` etc.)
###  
###  Notable packages that are NOT required: art, sbndcode, mrb.
###
###
################################################################################
cmake_minimum_required(VERSION 3.24)

# The presence of cetmodules appears to allow for some tricky UPS packages
# (e.g. `range` library) canvas depends on, which are not distributed
# with a package config script (neither pkgconfig or CMake flavor).
# Many other UPS packages appear they could work without it
# (hard to say since the dependency on `range` via `canvas` is quite hard)
find_package(cetmodules 3.20.00 REQUIRED)

#
# making `cmake` aware of Clang may require some care;
# we pick `gallery` qualifiers and if they contain `cN` then it's Clang;
# TODO check if this is still necessary with CMake 3.24
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
project(galleryAnalysis
  LANGUAGES CXX
  )

# this does not seem to affect the quality of the compilation,
# but it cuts a lot of the messages from CMake (and maybe it caches results?)
include(CetCMakeEnv)
cet_cmake_env()

# standard CMake voodoo to push the installation directory in the right direction
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}
    CACHE PATH "Install in source code dir" FORCE
  )
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


# this list can be filled last, from the target_link_libraries lists below
# most packages are looked for in CMAKE_PREFIX_PATH (UPS will populate it)
find_package(gallery         REQUIRED)
find_package(lardataalg      REQUIRED)
find_package(lardataobj      REQUIRED)
find_package(larcorealg      REQUIRED)
find_package(larcoreobj      REQUIRED)
find_package(sbndcode        REQUIRED)
find_package(fhiclcpp        REQUIRED)
find_package(messagefacility REQUIRED)
find_package(ROOT
  COMPONENTS Hist Tree RIO
  REQUIRED
  )


################################################################################
###  LIBRARY DIRECTIVES
################################################################################

#
# For this specific project, a library is not really needed: we could just link
# together all the source files and be gone with it.
# But hey, this is an example.
# 
# The INTERFACE keyword will make the symbols found in the libraries here
# implicitly visible to all targets depending on this one.
# Whether it is a good practice to use this feature this extensively is to be
# checked.
# 
# The list of linking targets was extracted from the headers included in the
# source code (`grep -h '#include' MCAssociations.* TrackAnalysis.* | sort -u`),
# although some of the correspondences are not straightforward (see comments).
#
add_library(AnalysisLib SHARED
  HitAnalysisAlg.cpp TrackAnalysis.cpp MCAssociations.cpp
  )
target_compile_features(AnalysisLib PUBLIC cxx_std_17)
target_link_libraries(AnalysisLib
  gallery::gallery
  lardataalg::headers        # for lardataalg/DetectorInfo/DetectorProperties.h
  larcorealg::Geometry       # for larcorealg/Geometry/GeometryCore.h
  lardataobj::AnalysisBase   # for lardataobj/AnalysisBase/BackTrackerMatchingData.h
  lardataobj::RecoBase       # for Track.h and Hit.h in lardataobj/RecoBase
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  canvas::canvas             # for canvas/Persistency/Common/FindMany.h
  ROOT::Hist
  ROOT::RIO
  ROOT::Core
  )
install(TARGETS AnalysisLib)


################################################################################
###  EXECUTABLE DIRECTIVES
################################################################################
add_executable(galleryAnalysis galleryAnalysis.cpp)
target_compile_features(galleryAnalysis PUBLIC cxx_std_17)
target_link_libraries(galleryAnalysis
  AnalysisLib
  gallery::gallery
  sbndcode::Geometry
  lardataalg::headers           # for headers in lardataalg/DetectorInfo
  lardataalg::DetectorInfo      # for lardataalg/DetectorInfo/DetectorPropertiesStandard.h
  larcorealg::Geometry          # for larcorealg/Geometry/GeometryCore.h
  lardataobj::RecoBase          # for lardataobj/RecoBase/Track.h
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
  ROOT::RIO
  )
install(TARGETS galleryAnalysis)