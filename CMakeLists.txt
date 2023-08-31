cmake_minimum_required(VERSION 3.18...3.22)

project(impa)

add_library(impalib INTERFACE)
target_include_directories(impalib INTERFACE include)

add_library(_wrapper SHARED src/impa/BCJR_wrapper_optimized.cpp)
target_link_libraries(_wrapper PRIVATE impalib)
install(TARGETS _wrapper LIBRARY DESTINATION src/impa)

add_subdirectory(test/external/cnpy)
add_executable(impalib_unit_tests test/src/impalib_unit_tests.cpp)
target_link_libraries(impalib_unit_tests PRIVATE impalib cnpy)
target_include_directories(impalib_unit_tests PRIVATE test/include test/external/cnpy)

add_executable(demo test/src/demo.cpp)
target_link_libraries(demo PRIVATE impalib)
