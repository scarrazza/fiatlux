noinst_PROGRAMS = generate
generate_SOURCES = generate.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -Wall -Wextra -fPIC -std=c++11 $(LHAPDF_CPPFLAGS)
AM_LDFLAGS = -L$(top_builddir)/src/.libs 
LIBS = -lfiatlux $(LHAPDF_LDFLAGS)
