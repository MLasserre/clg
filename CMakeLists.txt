cmake_minimum_required(VERSION 3.13)

project(
  CLG
  VERSION 0.1
  DESCRIPTION "CLG implementation using aGrUM and Eigen"
  LANGUAGES CXX)

set (CMAKE_CXX_STANDARD 14)


find_package(Eigen3 3.4 REQUIRED NO_MODULE)
if (Eigen3_FOUND)
    message (STATUS "Found Eigen3: ${Eigen3_DIR} (found version ${Eigen3_VERSION})")
    #set (TBB_LIBRARIES TBB::tbb)
else ()
    message(FATAL_ERROR "Please install Eigen3")
endif ()

find_package(aGrUM)
if (aGrUM_FOUND)
    message (STATUS "Found aGrUM: ${aGrUM_DIR} (found version ${aGrUM_VERSION})")
    include_directories(${AGRUM_INCLUDE_DIR})
    link_directories(${AGRUM_LIB_DIR})
else ()
    message(FATAL_ERROR "Please install aGrUM")
endif ()

add_subdirectory(src)
add_subdirectory(bin)
#add_subdirectory(include)


#if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    #set(CMAKE_CXX_EXTENSIONS OFF)
    #set_property(GLOBAL PROPERTY USE_FOLDERS ON)
    #include(CTest)
#endif()



#include_directories(src)

#add_subdirectory(src)

#if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR CLG_CMAKE_BUILD_TESTING)
    #AND BUILD_TESTING)
    #add_subdirectory(tests)
#endif()
