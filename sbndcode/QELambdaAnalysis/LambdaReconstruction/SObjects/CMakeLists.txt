if(DEFINED ENV{LARSOFT_VERSION})
    message(status " LARSOFT DEFINED SOBJECTS" )

    art_make_library(
        LIBRARY_NAME   sbndcode_SObjects
        
        SOURCE
            TPCSimpleHits.cpp
            TPCLinesHitDistanceUtils.cpp
            TPCSimpleClusters.cpp
            TPCSimpleLines.cpp
            TPCSimpleTriangles.cpp
            TPCLinesDistanceUtils.cpp
            TPCSimpleEvents.cpp
            TPCLinesPCA.cpp
            LambdaAnaTree.cpp

        LIBRARIES
            ROOT::Tree
            ROOT::Core
            ROOT::TMVA
      )

      install_headers()
      install_source()
endif()

message(status "Building " ${CMAKE_CURRENT_SOURCE_DIR})

add_library(sobjects_library
    TPCSimpleHits.cpp
    TPCSimpleHits.h
    TPCLinesHitDistanceUtils.cpp
    TPCLinesHitDistanceUtils.h
    TPCSimpleClusters.cpp
    TPCSimpleClusters.h
    TPCSimpleLines.cpp
    TPCSimpleLines.h
    TPCSimpleTriangles.cpp
    TPCSimpleTriangles.h
    TPCLinesDistanceUtils.cpp
    TPCLinesDistanceUtils.h
    TPCSimpleEvents.cpp
    TPCSimpleEvents.h
    TPCLinesPCA.cpp
    TPCLinesPCA.h
    LambdaAnaTree.cpp
    LambdaAnaTree.h
)

# Include directories for header files
target_include_directories(sobjects_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link the libraries
target_link_libraries(sobjects_library PUBLIC ROOT::Tree ROOT::Hist ROOT::Gpad)
#target_link_libraries(sobjects_library PUBLIC frans_library)