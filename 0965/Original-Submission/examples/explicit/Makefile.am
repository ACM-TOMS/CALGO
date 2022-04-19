AM_CXXFLAGS = --pedantic -Wall $(OPENMP_CXXFLAGS)
AM_CPPFLAGS = -I$(top_srcdir)/src

LOG_DRIVER = $(top_srcdir)/examples/explicit/run.sh

TESTS = explicit
check_PROGRAMS = explicit
explicit_SOURCES =  explicit.cpp
explicit_LDADD = $(top_builddir)/src/libridc.la 
