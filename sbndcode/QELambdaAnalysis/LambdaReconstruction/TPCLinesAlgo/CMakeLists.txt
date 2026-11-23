message(status "Building " ${CMAKE_CURRENT_SOURCE_DIR})

if(DEFINED ENV{LARSOFT_VERSION})
    message(status " LARSOFT DEFINED SOBJECTS" )

    art_make_library(
        LIBRARY_NAME   sbndcode_TPCLines
        
        SOURCE
            TPCLinesAlgo.cpp
            TPCSimpleCalo.cpp       
            TPCLinesDisplay.cpp            
            TPCLinesHough.cpp            
            TPCLinesTrackFinder.cpp            
            TPCLinesVertexFinder.cpp            
            TPCLinesDBSCAN.cpp            
            TPCLinesDirectionRecoUtils.cpp            
            TPCLinesParameters.cpp
        LIBRARIES
            ROOT::Tree
            ROOT::Core

            sbndcode_SObjects
            sbndcode_FRANS
      )

      install_headers()
      install_source()
endif()


add_library(tpclines_library
    TPCSimpleCalo.cpp
    TPCSimpleCalo.h
    TPCLinesAlgo.cpp
    TPCLinesAlgo.h
    TPCLinesDisplay.cpp
    TPCLinesDisplay.h
    TPCLinesHough.cpp
    TPCLinesHough.h
    TPCLinesTrackFinder.cpp
    TPCLinesTrackFinder.h
    TPCLinesVertexFinder.cpp
    TPCLinesVertexFinder.h
    TPCLinesDBSCAN.cpp
    TPCLinesDBSCAN.h
    TPCLinesDirectionRecoUtils.cpp
    TPCLinesDirectionRecoUtils.h
    TPCLinesParameters.cpp
    TPCLinesParameters.h
)


# Link math_library to string_library
target_link_libraries(tpclines_library PRIVATE ROOT::Hist ROOT::Graf ROOT::Gpad)
target_link_libraries(tpclines_library PRIVATE sobjects_library)
target_link_libraries(tpclines_library PRIVATE frans_library)


