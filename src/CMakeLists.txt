set(SOURCES
    canonical_form.cpp
    continuous_variable.cpp
    gaussian.cpp
    linear_gaussian.cpp
    scope.cpp
    variable_elimination.cpp
)

add_library(${PROJECT_NAME} ${SOURCES})


#add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include/CLG/
)

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE RELEASE)
endif()

target_link_libraries(${PROJECT_NAME} agrum Eigen3::Eigen)

