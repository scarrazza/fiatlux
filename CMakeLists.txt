cmake_minimum_required (VERSION 2.6)
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# the project name
project(fiatlux)

# activating some global properties for the project
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
enable_language(Fortran)

# export version to file
set(VERSION "\"0.0.2\"")
configure_file(
  "${PROJECT_SOURCE_DIR}/src/fiatlux/config.h.in"
  "${PROJECT_SOURCE_DIR}/src/fiatlux/config.h"
  )

set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")
configure_file(
  "${PROJECT_SOURCE_DIR}/scripts/fiatlux-config.in"
  "${PROJECT_SOURCE_DIR}/scripts/fiatlux-config"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/scripts/fiatlux.pc.in"
  "${PROJECT_SOURCE_DIR}/scripts/fiatlux.pc"
)

find_package(PkgConfig REQUIRED)
pkg_search_module(GSL REQUIRED gsl)
pkg_search_module(YAML REQUIRED yaml-cpp)

add_definitions(-DDATA_PATH="${prefix}/share/fiatlux")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${YAMLCPP_CFLAGS} ${GSL_CFLAGS} -std=c++11")
include_directories(src ${GSL_INCLUDE_DIRS} ${YAML_INCLUDE_DIRS})
add_library(fiatlux SHARED src/fiatlux.cc
                           src/elastic.cc
                           src/inelastic.cc
                           src/integrator.cc
                           src/msbar.cc
                           src/proton.cc
                           src/settings.cc
                           src/F2N.f
                           src/HERMES-ALLM-gd-fit-11.f
                         )

target_link_libraries(fiatlux ${GSL_LDFLAGS} ${YAML_LDFLAGS})

install(FILES ${PROJECT_SOURCE_DIR}/scripts/fiatlux-config
	DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(FILES ${PROJECT_SOURCE_DIR}/scripts/fiatlux.pc DESTINATION lib/pkgconfig)
install(DIRECTORY src/fiatlux DESTINATION include)
install(FILES ${PROJECT_SOURCE_DIR}/data/CrossSectionsOnly_SplinesWithVariableKnots.dat DESTINATION share/fiatlux)
install(FILES ${PROJECT_SOURCE_DIR}/data/CrossSectionsAndPolarized_SplinesWithVariableKnots.dat DESTINATION share/fiatlux)
install(TARGETS fiatlux DESTINATION lib)

option(ENABLE_EXAMPLES "Enable examples." OFF)
if(ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()
