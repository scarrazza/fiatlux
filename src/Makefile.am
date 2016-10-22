lib_LTLIBRARIES = libfiatlux.la
libfiatlux_la_SOURCES = fiatlux.cc \
			integrator.cc \
			elastic.cc \
			inelastic.cc \
			settings.cc \
			proton.cc

otherincludedir = $(includedir)/fiatlux
otherinclude_HEADERS = fiatlux/fiatlux.h \
		       fiatlux/integrator.h \
		       fiatlux/elastic.h \
		       fiatlux/inelastic.h \
		       fiatlux/settings.h \
		       fiatlux/proton.h

libfiatlux_la_LDFLAGS = -version-info 0:0:0 $(YAML_LDFLAGS)
libfiatlux_la_CPPFLAGS= -Wall -Wextra -fPIC $(YAML_CPPFLAGS)
