F77=$(FC)
FFLAGS=$(FCFLAGS)
FLIBS=$(FCLIBS)

LIBTIDES = $(top_builddir)/libTIDES/libTIDES.a

EXTRA_DIST = out_std_kepler_partials.dat out_std_lorenz_partials.dat

check_PROGRAMS = minc_kepler minc_lorenz

minc_kepler_SOURCES = test_minc_kepler.c minc_kepler.c minc_tides.c minc_tides.h
minc_lorenz_SOURCES = test_minc_lorenz.c minc_lorenz.c minc_tides.c minc_tides.h

if HAVE_FORTRAN_TESTS
  check_PROGRAMS += minf_kepler minf_lorenz
endif

minf_kepler_SOURCES = test_minf_kepler.f minf_kepler.f minf_tides.f
minf_lorenz_SOURCES = test_minf_lorenz.f minf_lorenz.f minf_tides.f

check_PROGRAMS += std_kepler std_kepler_partials

std_kepler_SOURCES = test_std_kepler.c std_kepler.c std_kepler.h

std_kepler_LDADD = $(LIBTIDES)
std_kepler_DEPENDENCIES = $(LIBTIDES)

std_kepler_partials_SOURCES = test_std_kepler_partials.c std_kepler_partials.c std_kepler_partials.h

std_kepler_partials_LDADD = $(LIBTIDES)
std_kepler_partials_DEPENDENCIES = $(LIBTIDES)

check_PROGRAMS += std_lorenz std_lorenz_partials

std_lorenz_SOURCES = test_std_lorenz.c std_lorenz.c std_lorenz.h

std_lorenz_LDADD = $(LIBTIDES)
std_lorenz_DEPENDENCIES = $(LIBTIDES)

std_lorenz_partials_SOURCES = test_std_lorenz_partials.c std_lorenz_partials.c std_lorenz_partials.h

std_lorenz_partials_LDADD = $(LIBTIDES)
std_lorenz_partials_DEPENDENCIES = $(LIBTIDES)

if HAVE_MP
  check_PROGRAMS += mpfr_kepler_14 mpfr_kepler_50 mpfr_lorenz_14 mpfr_lorenz_50
else
  LIBS += -lm
endif

mpfr_kepler_14_SOURCES = test_mpfr_kepler_14.c mpfr_kepler.c mpfr_kepler.h
mpfr_kepler_50_SOURCES = test_mpfr_kepler_50.c mpfr_kepler.c mpfr_kepler.h

mpfr_kepler_14_LDADD = $(LIBTIDES)
mpfr_kepler_14_DEPENDENCIES = $(LIBTIDES)
mpfr_kepler_50_LDADD = $(LIBTIDES)
mpfr_kepler_50_DEPENDENCIES = $(LIBTIDES)

mpfr_lorenz_14_SOURCES = test_mpfr_lorenz_14.c mpfr_lorenz.c mpfr_lorenz.h
mpfr_lorenz_50_SOURCES = test_mpfr_lorenz_50.c mpfr_lorenz.c mpfr_lorenz.h

mpfr_lorenz_14_LDADD = $(LIBTIDES)
mpfr_lorenz_14_DEPENDENCIES = $(LIBTIDES)
mpfr_lorenz_50_LDADD = $(LIBTIDES)
mpfr_lorenz_50_DEPENDENCIES = $(LIBTIDES)

TESTS = $(check_PROGRAMS)
