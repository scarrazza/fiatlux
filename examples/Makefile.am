noinst_PROGRAMS = generate
generate_SOURCES = generate.cc hoppet.f90

AM_CPPFLAGS = -I$(top_srcdir)/src -Wall -Wextra -fPIC -std=c++11 $(YAML_CPPFLAGS) $(APFEL_CPPFLAGS) $(LHAPDF_CPPFLAGS)
AM_FCFLAGS = $(HOPPET_FFFLAGS)
AM_LDFLAGS = -L$(top_builddir)/src/.libs 
LIBS = -lfiatlux $(HOPPET_LDFLAGS) $(YAML_LDFLAGS) $(APFEL_LDFLAGS) $(LHAPDF_LDFLAGS) -lgfortran
