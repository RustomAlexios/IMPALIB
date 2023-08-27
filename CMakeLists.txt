cmake_minimum_required(VERSION 3.18...3.22)

project(impa)

add_library(_wrapper SHARED src/impa/BCJR_wrapper_optimized.cpp)
target_include_directories(_wrapper PRIVATE include)
install(TARGETS _wrapper LIBRARY DESTINATION src/impa)
