cmake_minimum_required(VERSION 3.17)
project(sim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-no-pie")
set(CMAKE_C_FLAGS "-no-pie")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")


include_directories(algorithms)
add_subdirectory(algorithms lib)

add_executable(sim main.cpp)


target_link_libraries(sim boost_timer boost_random algorithms)