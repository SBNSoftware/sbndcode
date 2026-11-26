# find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Net)
find_package(Boost REQUIRED COMPONENTS system filesystem)


# using/not using larsoft
if(DEFINED ENV{LARSOFT_VERSION})
    message(status " LARSOFT DEFINED" )
    add_definitions(-DLAMBDAANA_LARSOFT=1)
else()
    message(status " LARSOFT NOT DEFINED" )
    add_definitions(-DLAMBDAANA_LARSOFT=0)
endif()



# Add subdirectories for libraries
add_subdirectory(EventHandle)
add_subdirectory(SObjects)
add_subdirectory(TPCLinesAlgo)
add_subdirectory(FRANS)
add_subdirectory(config)