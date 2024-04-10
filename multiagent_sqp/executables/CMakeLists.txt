add_executable(auv_sqp auv_sqp.cpp)

target_compile_features(auv_sqp PRIVATE cxx_std_20)

target_link_libraries(auv_sqp PRIVATE sqp_lib)
