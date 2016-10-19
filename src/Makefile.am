lib_LTLIBRARIES = libfiatlux.la
libfiatlux_la_SOURCES = fiatlux.cc \
			integrator.cc \
			elastic.cc \
			settings.cc

otherincludedir = $(includedir)/fiatlux
otherinclude_HEADERS = fiatlux/fiatlux.h \
		       fiatlux/integrator.h \
		       fiatlux/elastic.h \
		       fiatlux/settings.h

libfiatlux_la_LDFLAGS = -version-info 0:0:0 $(GSL_LDFLAGS)
libfiatlux_la_CPPFLAGS= -Wall -Wextra -fPIC $(GSL_CPPFLAGS)
