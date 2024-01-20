cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

# set C++ version requirements
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set variable for project name
set(PROJ_NAME
    stochastic-optimization-lib
)

# set libraries names (only libraries without templates)
set(LIBRARIES
)

# set project name and version
project(
  ${PROJ_NAME}
  VERSION 1.0
  LANGUAGES CXX)

# set include directory
include_directories(include src test)

# set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

# create all the executables
add_executable(test-saspso test/saspso.cpp)


# find openmp and eigen libraries
#list(APPEND CMAKE_PREFIX_PATH "")
find_package(OpenMP REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE
   HINTS $ENV{mkEigenPrefix}/share/eigen3/cmake
)
# If Eigen is not found:
#   Check if mk module eigen has been loaded
#   You may want to pass a cmake flag -DCMAKE_PREFIX_PATH=path/to/eigen/cmake/config

# link custom libraries, openmp and eigen to executables
if(OpenMP_CXX_FOUND)
    target_link_libraries(test-saspso PUBLIC  ${LIBRARIES} OpenMP::OpenMP_CXX Eigen3::Eigen)
else()
    target_link_libraries(test-saspso ${LIBRARIES})
endif()