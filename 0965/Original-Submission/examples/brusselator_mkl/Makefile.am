AM_CXXFLAGS = -Wall $(OPENMP_CXXFLAGS)
AM_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/src 


LOG_DRIVER = $(top_srcdir)/examples/brusselator_mkl/run.sh

TESTS = brusselator_mkl
check_PROGRAMS = brusselator_mkl
brusselator_mkl_SOURCES = brusselator.cpp  
brusselator_mkl_LDADD = $(top_builddir)/src/libridc.la -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
EXTRA_DIST = timing.sh
