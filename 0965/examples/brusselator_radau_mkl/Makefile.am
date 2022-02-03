AM_CXXFLAGS =  -Wall $(OPENMP_CXXFLAGS)
AM_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/src 


noinst_PROGRAMS = brusselator_radau
brusselator_radau_SOURCES = brusselator.cpp  radau.cpp
brusselator_radau_LDADD = -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
EXTRA_DIST = timing.sh ode.h
