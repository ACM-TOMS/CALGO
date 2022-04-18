AM_CXXFLAGS = -Wall $(OPENMP_CXXFLAGS)
AM_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/src 


LOG_DRIVER = $(top_srcdir)/examples/brusselator_gsl/run.sh

TESTS = brusselator_gsl
check_PROGRAMS = brusselator_gsl
brusselator_gsl_SOURCES = brusselator.cpp  
brusselator_gsl_LDADD = $(top_builddir)/src/libridc.la -lgsl -lm -lgslcblas
