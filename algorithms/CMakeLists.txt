cmake_minimum_required(VERSION 3.17)

project(algorithms)

file(GLOB_RECURSE LIB_FILES "*.h" "*.cpp")

add_library(algorithms ${LIB_FILES})

# set_target_properties(algorithms PROPERTIES LINKER_LANGUAGE CXX)
