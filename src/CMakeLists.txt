# Note that headers are optional, and do not affect add_library, but they will not
# show up in IDEs unless they are listed in add_library.

include(CMakePrintHelpers)

# Optionally glob, but only for CMake 3.12 or later:
file(GLOB_RECURSE HEADER_LIST CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../include/*.h")
cmake_print_variables(HEADER_LIST)

find_package(Eigen3 3.3 NO_MODULE)
find_package(SFML COMPONENTS graphics REQUIRED)
find_package(libInterpolate REQUIRED)


include(FetchContent)
FetchContent_Declare(
        tomlplusplus
        GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
        GIT_TAG        v3.3.0
)
FetchContent_MakeAvailable(tomlplusplus)

# TODO Private/Public headers.
# Make an automatic library - will be static or dynamic based on user setting
add_library(ma_interface ${HEADER_LIST} field_interface.cpp rrt_interface.cpp obstacle_interface.cpp sim_interface.cpp logger.cpp metrics.cpp bfs_planner.cpp visual.cpp linestrip.cpp)

#add_custom_target(copy_assets_ma
#        COMMAND ${CMAKE_COMMAND} -E copy_directory
#        ${CMAKE_CURRENT_SOURCE_DIR}/../resources/ $<TARGET_FILE_DIR:rrt_lib>)
#
#add_dependencies(ma_interace copy_assets_ma)

target_include_directories(ma_interface PUBLIC ../include ../resources PRIVATE ${tomlplusplus_SOURCE_DIR}/include)
target_link_libraries (ma_interface casadi sfml-graphics Eigen3::Eigen gmrf_lib rrt_lib sqp_lib random_walk_lib libInterpolate::Interpolate)


# All users of this library will need at least C++11
target_compile_features(ma_interface PUBLIC cxx_std_20)
# IDEs should put the headers in a nice place
source_group(
        TREE "${CMAKE_CURRENT_SOURCE_DIR}/../include"
        PREFIX "Header Files"
        FILES ${HEADER_LIST})

