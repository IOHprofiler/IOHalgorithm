#
cmake_minimum_required(VERSION 3.12) 


# set the project name and version
project(modularGA VERSION 1.1)

set(CMAKE_CXX_STANDARD 17)
set(EXTERNAL_DIR "${PROJECT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# set source file directorys and prepare for adding libraries.
file(GLOB_RECURSE HEADERS ./include *.hpp)
include_directories(GLOB_RECURSE)


add_library(modularGA INTERFACE)
add_library(modularGA::modularGA ALIAS modularGA)

# configure external package

include_directories(./IOH/include/ioh)
add_subdirectory(./IOH)
target_link_libraries(modularGA INTERFACE ioh)

include_directories(./IOH/external/fmt/include/)
target_link_libraries(modularGA INTERFACE fmt::fmt-header-only)

# Include external clutchlog lib
include_directories(./IOH/external/clutchlog)



target_include_directories(modularGA
INTERFACE
    # Interface paths are relative to the installation path, not the current directory.
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include> 
)

target_compile_features(modularGA
INTERFACE
    cxx_std_17
)


# Generate package config files.
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/cmake/modularGA-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/modularGA-config.cmake
INSTALL_DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/modularGA
) 

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/modularGA-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)


# Define installation rules.
include(GNUInstallDirs)

install(TARGETS modularGA fmt-header-only
    EXPORT modularGA-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
FILES
    ${CMAKE_CURRENT_BINARY_DIR}/modularGA-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/modularGA-config-version.cmake
DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/modularGA
)

install(EXPORT modularGA-targets
    FILE modularGA-targets.cmake
    NAMESPACE modularGA::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/modularGA
)
 

# Export targets. If an external project adds this project to its build tree
# via `add_subdirectory()`, then `find_package()` will work without
# requiring this project to be installed.
export(EXPORT modularGA-targets
    FILE ${CMAKE_CURRENT_BINARY_DIR}/modularGA-targets.cmake
    NAMESPACE modularGA::
)
export(PACKAGE modularGA)

# add the IOH library
# add the executable
add_subdirectory(run)