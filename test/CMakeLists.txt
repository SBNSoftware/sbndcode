# cmake driver file for testing from CET build tools
include(CetTest)

# Always enable asserts (for tests only)
cet_enable_asserts()

# test directories
add_subdirectory(Geometry)
add_subdirectory(LArSoftConfigurations)
add_subdirectory(JobConfigurations)
#add_subdirectory(CRT)
add_subdirectory(fcl)

# integration tests
add_subdirectory(ci)

