add_executable(manual manual.cpp)
add_executable(load_map load_map_test.cpp)
target_compile_features(manual PRIVATE cxx_std_20)
target_compile_features(load_map PRIVATE cxx_std_20)

target_link_libraries(manual PRIVATE gmrf_lib)
target_link_libraries(load_map PRIVATE gmrf_lib)
