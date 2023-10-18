cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(Plotting)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS RIO Net XMLIO)
find_program(LSB_RELEASE_EXEC lsb_release)
find_package(TBB REQUIRED)

include(${ROOT_USE_FILE})
string(CONCAT CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -Wall -Wextra -Wpedantic -Werror")
#grabs the distribution
execute_process(COMMAND ${LSB_RELEASE_EXEC} -is
  OUTPUT_VARIABLE LSB_RELEASE_ID_SHORT
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

string(FIND "${LSB_RELEASE_ID_SHORT}" "Scientific" _cmp)
if (_cmp GREATER -1)
  set(Boost_NO_BOOST_CMAKE ON)
else()
  string(CONCAT CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -std=c++1y")
endif()

find_package(Boost REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIRS})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/IO)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Analyzer)

