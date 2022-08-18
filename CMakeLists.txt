# CMakeLists.txt
cmake_minimum_required(VERSION 3.0)
project(test)
find_package(Eigen3 CONFIG REQUIRED)
add_executable(NNA NNA.cpp)
target_link_libraries(NNA PRIVATE Eigen3::Eigen)
