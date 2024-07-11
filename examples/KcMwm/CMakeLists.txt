cmake_minimum_required(VERSION 3.18...3.22)

project(impalib-demo VERSION "0.1.0" LANGUAGES CXX)

list(APPEND CMAKE_PREFIX_PATH ${IMPALIB_CMAKE_PACKAGE_INSTALL_SUBDIR})

find_package(impalib REQUIRED)
add_executable(demo demo.cpp)
target_link_libraries(demo PRIVATE impalib::impalib)
