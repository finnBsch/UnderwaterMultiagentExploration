cmake_minimum_required(VERSION 3.21)

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fno-stack-protector")

# Project name and a few useful settings. Other commands can pick up the results
project(
        i_rrtstar
        VERSION 0.1
        DESCRIPTION "i_rrtstar"
        LANGUAGES CXX)


# FetchContent added in CMake 3.11, downloads during the configure step
# FetchContent_MakeAvailable was added in CMake 3.14; simpler usage
include(FetchContent)

# Accumulator library
# This is header only, so could be replaced with git submodules or FetchContent
find_package(Boost REQUIRED)
# Adds Boost::boost


# The compiled library code is here
add_subdirectory(src)

# The executable code is here
add_subdirectory(executables)
