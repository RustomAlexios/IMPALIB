cmake_minimum_required(VERSION 3.18...3.22)

project(impa)

add_library(impalib INTERFACE)
target_include_directories(impalib INTERFACE include)

add_library(_wrapper SHARED src/impa/BCJR_wrapper_optimized.cpp)
target_link_libraries(_wrapper PRIVATE impalib)
install(TARGETS _wrapper LIBRARY DESTINATION src/impa)

