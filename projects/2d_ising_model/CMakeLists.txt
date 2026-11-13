cmake_minimum_required(VERSION 3.18)

project(IsingModel LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED NO_MODULE)

add_executable(IsingModel main.cpp)

target_link_libraries(IsingModel PRIVATE Eigen3::Eigen)