cmake_minimum_required(VERSION 3.10)
project(robot_camera_calibration)

set(CMAKE_CXX_STANDARD 11)

find_package(GSL REQUIRED)

add_executable(robot_camera_calibration main.cpp utils.cpp utils.h Global_parameters.h)

target_link_libraries(robot_camera_calibration GSL::gsl GSL::gslcblas)