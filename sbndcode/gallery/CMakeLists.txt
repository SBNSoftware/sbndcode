add_subdirectory(helpers)
add_subdirectory(python)

foreach(ExampleDir IN ITEMS galleryAnalysis )
  
  # take *all* the files from the listed directory
  file(GLOB_RECURSE GalleryExamples RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${ExampleDir}/* )
  
  # install will reproduce the intermediate subdirectories from the source list 
  # into the destination directory
  install(FILES ${GalleryExamples} DESTINATION "examples/gallery/${ExampleDir}")
  
endforeach()
