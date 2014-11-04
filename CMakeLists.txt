cmake_minimum_required(VERSION 2.8.7)
project(gamgee)

include(ExternalProject)

## This is the right way to do this, but only supported in cmake 2.8.12
#add_compile_options("-std=c++1y")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

# Dependency: Boost Unit Test Framework (find in the system)
find_package(Boost 1.55 COMPONENTS unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# enable installing dependencies
option(INSTALL_DEPENDENCIES
       "Install project dependencies"
       OFF)

# Dependency: htslib (download and build)
include("contrib/htslib.cmake")

include_directories(gamgee)

add_subdirectory(gamgee)
add_subdirectory(test)

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )
