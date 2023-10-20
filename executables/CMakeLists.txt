include(CMakePrintHelpers)
add_executable(full_interface full_interface.cpp)

target_compile_features(full_interface PRIVATE cxx_std_20)


find_package(OpenMP)

target_link_libraries(full_interface PRIVATE rrt_lib sqp_lib ma_interface OpenMP::OpenMP_CXX)

