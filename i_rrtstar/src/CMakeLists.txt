# Note that headers are optional, and do not affect add_library, but they will not
# show up in IDEs unless they are listed in add_library.

include(CMakePrintHelpers)

set(RES_LOC ${CMAKE_BINARY_DIR}/../resources_build)


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../include/res_loc.h.in ${CMAKE_CURRENT_SOURCE_DIR}/../include/res_loc.h @ONLY)

# Optionally glob, but only for CMake 3.12 or later:
file(GLOB_RECURSE HEADER_LIST CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../include/*.h")
cmake_print_variables(HEADER_LIST)

find_package(Eigen3 3.3 NO_MODULE)
find_package(SFML COMPONENTS graphics REQUIRED)


include(FetchContent)
FetchContent_Declare(
        tomlplusplus
        GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
        GIT_TAG        v3.3.0
)
FetchContent_MakeAvailable(tomlplusplus)

# TODO Private/Public headers.
# Make an automatic library - will be static or dynamic based on user setting
add_library(rrt_lib  ${HEADER_LIST} util.cpp viz.cpp rrt.cpp rrt_configs.cpp)

add_custom_target(copy_assets_rrt
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/../resources/ ${RES_LOC})

add_dependencies(rrt_lib copy_assets_rrt)

target_include_directories(rrt_lib PUBLIC ../include ../resources PRIVATE ${tomlplusplus_SOURCE_DIR}/include)
target_link_libraries (rrt_lib casadi sfml-graphics Eigen3::Eigen)


# All users of this library will need at least C++11
target_compile_features(rrt_lib PUBLIC cxx_std_20)
# IDEs should put the headers in a nice place
source_group(
        TREE "${CMAKE_CURRENT_SOURCE_DIR}/../include"
        PREFIX "Header Files"
        FILES ${HEADER_LIST})

