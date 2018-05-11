#
# this code applies to all the directories listed in the first line:
# - processes all the ".in" files replacing the @XXX@ placeholders
#   - the new files are added to the list of things to install
# - installs all the files, where the ".in" files are replaced by the processed ones
#

if(NOT CMAKE_CXX_STANDARD)
  # CMAKE_CXX_STANDARD is not set: figure it out
  string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" FlagsVarName)
  string(REGEX REPLACE "(^|.* )-std=c\\+\\+([^ ]*)( .*|$)" "\\2" MyStandardsAreBetter "${${FlagsVarName}}")
else()
  set(MyStandardsAreBetter "$CMAKE_CXX_STANDARD")
endif()

foreach(ExampleDir IN ITEMS galleryAnalysis )
  
  # take *all* the files from the listed directory
  file(GLOB_RECURSE GalleryExamples RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${ExampleDir}/* )
  file(GLOB_RECURSE GalleryExamplesInput RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${ExampleDir}/*.in )
  
  foreach(ConfigInputFile IN LISTS GalleryExamplesInput )
    
    string(REGEX REPLACE ".in$" "" ConfigOutputFile ${ConfigInputFile})
    string(PREPEND ConfigOutputFile "${CMAKE_CURRENT_BINARY_DIR}/" )
    configure_file("${ConfigInputFile}" "${ConfigOutputFile}" @ONLY)
    list(REMOVE_ITEM GalleryExamples "${ConfigInputFile}")
    list(APPEND GalleryExamples "${ConfigOutputFile}")
    
  endforeach()
  
  # cet_install_files() will reproduce the intermediate subdirectories from the
  # source list into the destination directory
  cet_install_files(
    LIST "${GalleryExamples}"
    DIRNAME "examples/gallery/${ExampleDir}"
    )
  
endforeach()

#
# back to the routine:
#
add_subdirectory(helpers)
