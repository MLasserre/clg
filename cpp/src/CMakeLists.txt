cmake_minimum_required(VERSION 3.0)
project(CLG)

set (CMAKE_CXX_STANDARD 14)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

set(AGRUM_INSTALLATION_DIRECTORY "/home/lasserre/.venv/otagrum")
set(aGrUM_DIR "${AGRUM_INSTALLATION_DIRECTORY}/home/lasserre/.venv/otagrum/lib/cmake/aGrUM/")

find_package(aGrUM)
if (aGrUM_FOUND)
  include_directories(${AGRUM_INCLUDE_DIR})
  link_directories(${AGRUM_LIB_DIR})
else (aGrUM_FOUND)
  message(FATAL_ERROR "Please install aGrUM")
endif (aGrUM_FOUND)

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE RELEASE)
endif()

file(GLOB CLG_SOURCE ${CLG_SOURCE_DIR}/*.cpp)
file(GLOB CLG_INCLUDE ${CLG_SOURCE_DIR}/*.h)
add_executable (lauritzen_1992_example ${CLG_SOURCE})
if (${CMAKE_BUILD_TYPE}  STREQUAL "RELEASE") # release : act install release
  target_link_libraries(lauritzen_1992_example agrum Eigen3::Eigen)
else() # debug : act install debug
  target_link_libraries(lauritzen_1992_example agrum-dbg)
endif()
