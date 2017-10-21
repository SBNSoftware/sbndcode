add_subdirectory(helpers)

foreach(ExampleDir IN ITEMS galleryAnalysis )
  
  # take *all* the files from the listed directory
  file(GLOB_RECURSE GalleryExamples RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${ExampleDir}/* )
  
  # cet_install_files() will reproduce the intermediate subdirectories from the
  # source list into the destination directory
  cet_install_files(
    LIST ${GalleryExamples}
    DIRNAME examples/gallery/${ExampleDir}
    )
  
endforeach()
