cmake_minimum_required(VERSION 3.18...3.22)

project(impalib-demo VERSION "0.1.0" LANGUAGES CXX)

find_package(impalib REQUIRED)
add_executable(demo demo.cpp)
target_link_libraries(demo PRIVATE impalib::impalib)
