
add_executable(auv auv_control.cpp ../include/info_field.h ../include/information_cost.h ../include/path_tracking.h ../include/path.h ../include/com_path_tracking.h ../include/corridor_constraint.h )
add_executable(auv_low_level auv_control_low_level.cpp)
add_executable(auv_sqp auv_sqp.cpp)


target_compile_features(auv PRIVATE cxx_std_20)
target_compile_features(auv_low_level PRIVATE cxx_std_20)
target_compile_features(auv_sqp PRIVATE cxx_std_20)

target_link_libraries(auv PRIVATE sqp_lib)
target_link_libraries(auv_low_level PRIVATE sqp_lib)
target_link_libraries(auv_sqp PRIVATE sqp_lib)
