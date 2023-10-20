cmake_minimum_required(VERSION 3.21)

project(
        MASTERMAIN
        VERSION 0.1
        DESCRIPTION "MASTERMAIN"
        LANGUAGES CXX)


include(FetchContent)

find_package(Boost REQUIRED)



add_subdirectory(multiagent_sqp/executables)
#add_subdirectory(i_rrtstar/executables)
add_subdirectory(executables)
add_subdirectory(gmrf/executables)
add_subdirectory(i_rrtstar/executables)

add_subdirectory(src)
add_subdirectory(multiagent_sqp/src)
add_subdirectory(i_rrtstar/src)
add_subdirectory(gmrf/src)
add_subdirectory(random_walks/src)


