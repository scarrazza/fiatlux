cmake_minimum_required (VERSION 3.15..3.22)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# the project name
project(fiatlux VERSION "0.1.0")

# activating some global properties for the project
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC")
enable_language(Fortran)

set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")
configure_file(
    "${PROJECT_SOURCE_DIR}/src/fiatlux/config.h.in"
    "${PROJECT_SOURCE_DIR}/src/fiatlux/config.h"
    )

add_definitions(-DDATA_PATH="${prefix}/share/fiatlux")
include_directories(src)
add_subdirectory(external/yaml-cpp)

option(PYTHON_ONLY "Compiles only for python." OFF)
if(PYTHON_ONLY)
  if(SKBUILD)
    execute_process(
      COMMAND "${PYTHON_EXECUTABLE}" -c
              "import pybind11; print(pybind11.get_cmake_dir())"
      OUTPUT_VARIABLE _tmp_dir
      OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
    list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
  endif()
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_core MODULE src/pywrapper.cc
                           src/fiatlux.cc
                           src/elastic.cc
                           src/inelastic.cc
                           src/integrator.cc
                           src/msbar.cc
                           src/proton.cc
                           src/settings.cc
                           src/F2N.f
                           src/HERMES-ALLM-gd-fit-11.f)
  target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
  target_link_libraries(_core PRIVATE yaml-cpp)
  install(FILES ${PROJECT_SOURCE_DIR}/data/CrossSectionsOnly_SplinesWithVariableKnots.dat DESTINATION share/fiatlux)
  install(FILES ${PROJECT_SOURCE_DIR}/data/CrossSectionsAndPolarized_SplinesWithVariableKnots.dat DESTINATION share/fiatlux)
  install(TARGETS _core DESTINATION .)
else()
  configure_file(
    "${PROJECT_SOURCE_DIR}/scripts/fiatlux-config.in"
    "${PROJECT_SOURCE_DIR}/scripts/fiatlux-config"
  )
  configure_file(
    "${PROJECT_SOURCE_DIR}/scripts/fiatlux.pc.in"
    "${PROJECT_SOURCE_DIR}/scripts/fiatlux.pc"
  )
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
  target_link_libraries(fiatlux yaml-cpp)

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
endif()
