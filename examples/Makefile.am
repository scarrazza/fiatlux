noinst_PROGRAMS = generate
generate_SOURCES = generate.cc \
		   hoppet-struct-fns/xc3ns2p.f \
		   hoppet-struct-fns/xc2ns2p.f \
		   hoppet-struct-fns/xc2sg2p.f \
		   hoppet-struct-fns/xclns2p.f \
		   hoppet-struct-fns/xclsg2p.f \
		   hoppet-struct-fns/xc3ns3p.f \
		   hoppet-struct-fns/xcdiff3p.f \
		   hoppet-struct-fns/xc2ns3p.f \
		   hoppet-struct-fns/xc2sg3p.f \
		   hoppet-struct-fns/xclns3p.f \
		   hoppet-struct-fns/xclsg3p.f \
		   hoppet-struct-fns/io_utils.f90 \
		   hoppet-struct-fns/lcl_dec.f90 \
		   hoppet-struct-fns/parameters.f90 \
		   hoppet-struct-fns/helpers.f90 \
		   hoppet-struct-fns/hoppet_coefficient_functions_N3LO.f90 \
		   hoppet-struct-fns/hoppet_structure_functions_N3LO.f90 \
		   hoppet.f90

AM_CPPFLAGS = -I$(top_srcdir)/src -Wall -Wextra -fPIC -std=c++11 $(YAML_CPPFLAGS) $(APFEL_CPPFLAGS) $(LHAPDF_CPPFLAGS)
AM_FCFLAGS = $(HOPPET_FFFLAGS) -I./
AM_LDFLAGS = -L$(top_builddir)/src/.libs 
LIBS = -lfiatlux $(HOPPET_LDFLAGS) $(YAML_LDFLAGS) $(APFEL_LDFLAGS) $(LHAPDF_LDFLAGS) -lgfortran 
