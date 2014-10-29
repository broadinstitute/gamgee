cmake_minimum_required(VERSION 2.8.4)
project(gamgee)

include(ExternalProject)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

find_package(Boost 1.56 COMPONENTS unit_test_framework REQUIRED)

include_directories(gamgee)
include_directories(gamgee/utils)
include_directories(lib)
include_directories(lib/htslib)

add_subdirectory(gamgee)
add_subdirectory(test)
