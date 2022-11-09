# ======================================================================
#  larsoft main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsoft/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


cmake_minimum_required(VERSION 3.20 FATAL_ERROR) ### MIGRATE-CHANGED (migrate-3.20.00) - VERSION 3.19 -> 3.20

set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 09.62.00)
find_package(cetmodules REQUIRED) ### MIGRATE-ADDED (migrate-3.20.00) - find_package(cetmodules) must precede project()
project(sbndcode LANGUAGES CXX)

# for CI:
# project(sbndcode VERSION 09.62.00)

message(STATUS "\n")
message(STATUS "=================================   ${PROJECT_NAME}   =================================")
message(STATUS "")


include(CetCMakeEnv)
cet_cmake_env()

# The specification of the C++ standard is required because CET macros
# set it to C++98 by default, and then override it depending on the
# compiler version, as determined by the compiler qualifier.
# Removing that qualifier, we need to supply that information manually,
# and maintain it too.
cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)
##add_definitions(-DGENIE_PRE_R3)

cet_report_compiler_flags() ### MIGRATE-ACTION-RECOMMENDED (migrate-3.20.00) - add args: REPORT_THRESHOLD VERBOSE

# save the repository tag
unset(GIT_REPOSITORY_TAG)
execute_process(
  COMMAND git -C "${CMAKE_CURRENT_SOURCE_DIR}" describe --dirty
  OUTPUT_VARIABLE GIT_REPOSITORY_TAG
  )
file(WRITE RepositoryTag "${GIT_REPOSITORY_TAG}")
message(STATUS "${PROJECT_NAME} repository tag: ${GIT_REPOSITORY_TAG}")
install_source(LIST RepositoryTag)

# these are minimum required versions, not the actual product versions
find_package( sbncode REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(sbncode ...) -> find_package(sbncode REQUIRED ...)
find_package( sbnobj REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(sbnobj ...) -> find_package(sbnobj REQUIRED ...)
find_package( larpandora REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(larpandora ...) -> find_package(larpandora REQUIRED ...)
find_package( larsim REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(larsim ...) -> find_package(larsim REQUIRED ...)
find_package( lardata REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(lardata ...) -> find_package(lardata REQUIRED ...)
find_package( lardataobj REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(lardataobj ...) -> find_package(lardataobj REQUIRED ...)
find_package( lardataalg REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(lardataalg ...) -> find_package(lardataalg REQUIRED ...)
find_package( larevt REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(larevt ...) -> find_package(larevt REQUIRED ...)
find_package( larcore REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(larcore ...) -> find_package(larcore REQUIRED ...)
find_package( larcorealg REQUIRED )
find_package( larcoreobj REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(larcoreobj ...) -> find_package(larcoreobj REQUIRED ...)
find_package( larreco REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(larreco ...) -> find_package(larreco REQUIRED ...)
find_package( larana REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(larana ...) -> find_package(larana REQUIRED ...)
find_package( nutools REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(nutools ...) -> find_package(nutools REQUIRED ...)
find_package( nusimdata REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(nusimdata ...) -> find_package(nusimdata REQUIRED ...)
find_package( nug4 REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(nug4 ...) -> find_package(nug4 REQUIRED ...)
find_package( nugen REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(nugen ...) -> find_package(nugen REQUIRED ...)
find_package( nuevdb REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(nuevdb ...) -> find_package(nuevdb REQUIRED ...)
find_package( nurandom REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(nurandom ...) -> find_package(nurandom REQUIRED ...)
find_package( dk2nudata )
find_package( dk2nugenie )
find_package( art REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(art ...) -> find_package(art REQUIRED ...)
find_package( art_root_io REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(art_root_io ...) -> find_package(art_root_io REQUIRED ...)
find_package( PostgreSQL REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(postgresql ...) -> find_package(postgresql REQUIRED ...)
find_package( hep_concurrency REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(hep_concurrency ...) -> find_package(hep_concurrency REQUIRED ...)
find_package( Eigen3 REQUIRED ) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(eigen ...) -> find_package(eigen REQUIRED ...)
find_package(Geant4 REQUIRED) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(...) -> find_package(REQUIRED ...)
find_package(Boost REQUIRED) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(...) -> find_package(REQUIRED ...)
find_package(ROOT REQUIRED) ### MIGRATE-CHANGED (migrate-3.20.00) - find_ups_product(...) -> find_package(REQUIRED ...)
find_package(artdaq_core REQUIRED)
find_package(sbndaq_artdaq_core REQUIRED)

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(sbndcode)
# subdirectory for test code
add_subdirectory(test)
# ups - table and config files
# subdirectory for scripts
add_subdirectory(bin)

# packaging utility

cet_cmake_config() ### MIGRATE-ADDED (migrate-3.20.00) - required cmd
