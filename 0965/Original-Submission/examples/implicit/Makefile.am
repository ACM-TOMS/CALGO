AM_CXXFLAGS = --pedantic -Wall $(OPENMP_CXXFLAGS)
AM_CPPFLAGS = -I$(top_srcdir)/src

LOG_DRIVER = $(top_srcdir)/examples/implicit/run.sh

TESTS = implicit
check_PROGRAMS = implicit
implicit_SOURCES = implicit.cpp  
implicit_LDADD = $(top_builddir)/src/libridc.la 
