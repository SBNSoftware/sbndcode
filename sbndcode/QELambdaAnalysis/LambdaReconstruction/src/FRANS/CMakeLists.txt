message(status "Building " ${CMAKE_CURRENT_SOURCE_DIR})


if(DEFINED ENV{LARSOFT_VERSION})
    message(status " LARSOFT DEFINED FRANS" )

    art_make_library( 
        LIBRARY_NAME   sbndcode_FRANS
        
        SOURCE
          FRANSObj.cpp
          FRANSTree.cpp
          ChargeDensity.cpp
          LeastSquares.h
          VertexView.h
        
        LIBRARIES
        ROOT::Tree
        ROOT::Core
        ROOT::TMVA

        sbndcode_SObjects
      
      )

      art_dictionary(DICTIONARY_LIBRARIES sbndcode_FRANS)

      install_headers()
      install_source()
endif()

add_library(frans_library
  LeastSquares.h
  FRANSObj.cpp
  FRANSObj.h
  ChargeDensity.cpp
  ChargeDensity.h
  ChargeDensityPset.h
  FRANSTree.cpp
  FRANSTree.h
  VertexView.h
)

# Include directories for header files
target_include_directories(frans_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link the libraries
target_link_libraries(frans_library PRIVATE sobjects_library)
target_link_libraries(frans_library PRIVATE ROOT::Hist ROOT::Gpad)
target_link_libraries(frans_library PUBLIC ROOT::TMVA)
