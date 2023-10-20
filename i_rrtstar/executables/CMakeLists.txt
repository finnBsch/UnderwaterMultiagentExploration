add_executable(rrt_moving rrt_moving.cpp)
add_executable(rrt_static rrt_static.cpp)

target_compile_features(rrt_moving PRIVATE cxx_std_20)
target_compile_features(rrt_static PRIVATE cxx_std_20)

target_link_libraries(rrt_moving PRIVATE rrt_lib)
target_link_libraries(rrt_static PRIVATE rrt_lib)
