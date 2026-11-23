message(status "Building " ${CMAKE_CURRENT_SOURCE_DIR})

if(DEFINED ENV{LARSOFT_VERSION})
    message(status " LARSOFT DEFINED SOBJECTS" )

    art_make_library(
        LIBRARY_NAME   sbndcode_SEventHandle
        
        SOURCE
            CommandLineParser.cpp
            SEventHandle.cpp
            STPCAnalyzerTreeReader.cpp
        LIBRARIES
            ROOT::Tree
            ROOT::Core

            sbndcode_SObjects
            sbnobj::Common_CRT
      )

      install_headers()
      install_source()
endif()

add_library(seventhandle_library
    CommandLineParser.cpp
    CommandLineParser.h
    SEventHandle.cpp
    SEventHandle.h
    STPCAnalyzerTreeReader.cpp
    STPCAnalyzerTreeReader.h
    SParserConfig.h
)

# Include directories for header files
target_include_directories(seventhandle_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# Link the libraries
target_link_libraries(seventhandle_library PUBLIC frans_library tpclines_library sobjects_library ROOT::Tree  ROOT::Tree Boost::filesystem)