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
find_package(casadi)

include(FetchContent)
FetchContent_Declare(
        tomlplusplus
        GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
        GIT_TAG        v3.3.0
)
find_package(OpenMP)

FetchContent_MakeAvailable(tomlplusplus)

# TODO Private/Public headers.
# Make an automatic library - will be static or dynamic based on user setting
add_library(sqp_lib configs.cpp connectivity_constraint.cpp bicycle_dynamics.cpp sqp_multiagent.cpp sqp_module.cpp ${HEADER_LIST}  collision_constraint.cpp information_cost.cpp path.cpp com_path_tracking.cpp corridor_constraint.cpp info_field.cpp ../include/hippocampus_dynamics.h hippocampus_dynamics.cpp bluerov_dynamics.cpp distribution_cost.cpp opti_wrapper.cpp nearest_obstacle_avoidance.cpp rrt_util.cpp)
#add_library(solve_compiled SHARED solver/solve_compiled.cpp solver/solve_compiled.h)

add_custom_target(copy_assets_sqp
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/../resources/ ${RES_LOC})

add_dependencies(sqp_lib copy_assets_sqp)

target_include_directories(sqp_lib PUBLIC ../include PRIVATE ${tomlplusplus_SOURCE_DIR}/include)
#target_include_directories(solve_compiled PUBLIC solver)
target_link_libraries (sqp_lib casadi sfml-graphics Eigen3::Eigen osqp  OpenMP::OpenMP_CXX)
#target_link_libraries(solve_compiled osqp)


# All users of this library will need at least C++11
target_compile_features(sqp_lib PUBLIC cxx_std_20)
# IDEs should put the headers in a nice place
source_group(
        TREE "${CMAKE_CURRENT_SOURCE_DIR}/../include"
        PREFIX "Header Files"
        FILES ${HEADER_LIST})

