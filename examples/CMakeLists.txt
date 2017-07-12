# LHAPDF
find_program(LHAPDF_CONFIG lhapdf-config REQUIRED)
if (LHAPDF_CONFIG)
  exec_program(${LHAPDF_CONFIG}
    ARGS --cflags
    OUTPUT_VARIABLE LHAPDF_CXX_FLAGS
  )
  set(LHAPDF_CXX_FLAGS ${LHAPDF_CXX_FLAGS} CACHE STRING INTERNAL)
  exec_program(${LHAPDF_CONFIG}
    ARGS --libs
    OUTPUT_VARIABLE LHAPDF_LIBRARIES
  )
  set(LHAPDF_LIBRARIES ${LHAPDF_LIBRARIES} CACHE STRING INTERNAL)
endif(LHAPDF_CONFIG)

# APFEL
find_program(APFEL_CONFIG apfel-config REQUIRED)
if (APFEL_CONFIG)
  exec_program(${APFEL_CONFIG}
    ARGS --cppflags
    OUTPUT_VARIABLE APFEL_CXX_FLAGS
  )
  set(APFEL_CXX_FLAGS ${APFEL_CXX_FLAGS} CACHE STRING INTERNAL)
  exec_program(${APFEL_CONFIG}
    ARGS --ldflags
    OUTPUT_VARIABLE APFEL_LIBRARIES
  )
  set(APFEL_LIBRARIES ${APFEL_LIBRARIES} CACHE STRING INTERNAL)
endif(APFEL_CONFIG)

# HOPPET
find_program(HOPPET_CONFIG hoppet-config REQUIRED)
if (HOPPET_CONFIG)
  exec_program(${HOPPET_CONFIG}
    ARGS --cxxflags
    OUTPUT_VARIABLE HOPPET_CXX_FLAGS
  )
  set(HOPPET_CXX_FLAGS ${HOPPET_CXX_FLAGS} CACHE STRING INTERNAL)
  exec_program(${HOPPET_CONFIG}
    ARGS --ldflags
    OUTPUT_VARIABLE HOPPET_LIBRARIES
  )
  set(HOPPET_LIBRARIES ${HOPPET_LIBRARIES} CACHE STRING INTERNAL)
endif(HOPPET_CONFIG)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LHAPDF_CXX_FLAGS} ${APFEL_CXX_FLAGS} ${HOPPET_CXX_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${HOPPET_CXX_FLAGS}/hoppet")

add_library(hoppet_bridge STATIC
  hoppet-struct-fns/xc3ns2p.f
  hoppet-struct-fns/xc2ns2p.f
  hoppet-struct-fns/xc2sg2p.f
  hoppet-struct-fns/xclns2p.f
  hoppet-struct-fns/xclsg2p.f
  hoppet-struct-fns/xc3ns3p.f
  hoppet-struct-fns/xcdiff3p.f
  hoppet-struct-fns/xc2ns3p.f
  hoppet-struct-fns/xc2sg3p.f
  hoppet-struct-fns/xclns3p.f
  hoppet-struct-fns/xclsg3p.f
  hoppet-struct-fns/io_utils.f90
  hoppet-struct-fns/lcl_dec.f90
  hoppet-struct-fns/parameters.f90
  hoppet-struct-fns/helpers.f90
  hoppet-struct-fns/hoppet_coefficient_functions_N3LO.f90
  hoppet-struct-fns/hoppet_structure_functions_N3LO.f90
  hoppet.f90
)

add_executable(comparealpha comparealpha.cc)
target_link_libraries(comparealpha fiatlux hoppet_bridge ${APFEL_LIBRARIES} ${HOPPET_LIBRARIES} ${LHAPDF_LIBRARIES})

add_executable(comparesf comparesf.cc)
target_link_libraries(comparesf fiatlux hoppet_bridge ${APFEL_LIBRARIES} ${HOPPET_LIBRARIES} ${LHAPDF_LIBRARIES})

add_executable(generate generate.cc)
target_link_libraries(generate fiatlux hoppet_bridge ${APFEL_LIBRARIES} ${HOPPET_LIBRARIES} ${LHAPDF_LIBRARIES})
