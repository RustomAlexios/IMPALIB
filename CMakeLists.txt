cmake_minimum_required(VERSION 3.18...3.22)

project(impalib VERSION "0.1.0" LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)

include(GNUInstallDirs)

add_library(impalib INTERFACE)
target_include_directories(
        impalib INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS impalib EXPORT impalibTargets)

set(IMPALIB_CMAKE_PACKAGE_INSTALL_SUBDIR "share/impalib/cmake")

install(
        EXPORT impalibTargets
        NAMESPACE impalib::
        DESTINATION ${IMPALIB_CMAKE_PACKAGE_INSTALL_SUBDIR})

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
        impalibConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMinorVersion)

configure_package_config_file(
        "${PROJECT_SOURCE_DIR}/cmake/impalibConfig.cmake.in" impalibConfig.cmake
        INSTALL_DESTINATION ${IMPALIB_CMAKE_PACKAGE_INSTALL_SUBDIR})

install(FILES "${PROJECT_BINARY_DIR}/impalibConfig.cmake"
        "${PROJECT_BINARY_DIR}/impalibConfigVersion.cmake"
        DESTINATION ${IMPALIB_CMAKE_PACKAGE_INSTALL_SUBDIR})

if(DEFINED SKBUILD)
        add_library(_wrapper SHARED src/impa/Wrapper.cpp)
        target_link_libraries(_wrapper PRIVATE impalib)
        install(TARGETS _wrapper LIBRARY DESTINATION src/impa)

else()
        add_subdirectory(test/external/cnpy)
        add_executable(impalib_unit_tests_kc_mwm test/src/impalib_unit_tests_kc_mwm.cpp)
        target_link_libraries(impalib_unit_tests_kc_mwm PRIVATE impalib cnpy)
        target_include_directories(impalib_unit_tests_kc_mwm PRIVATE test/include test/external/cnpy)

        add_executable(impalib_unit_tests_tsp test/src/impalib_unit_tests_tsp.cpp)
        target_link_libraries(impalib_unit_tests_tsp PRIVATE impalib cnpy)
        target_include_directories(impalib_unit_tests_tsp PRIVATE test/include test/external/cnpy)

        add_executable(impalib_unit_tests_ksat test/src/impalib_unit_tests_ksat.cpp)
        target_link_libraries(impalib_unit_tests_ksat PRIVATE impalib cnpy)
        target_include_directories(impalib_unit_tests_ksat PRIVATE test/include test/external/cnpy)

        file(COPY test DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
        file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
