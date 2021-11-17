noinst_LTLIBRARIES = libgslwavelet.la 

pkginclude_HEADERS = gsl_wavelet.h gsl_wavelet2d.h

INCLUDES= -I$(top_builddir)

libgslwavelet_la_SOURCES = dwt.c wavelet.c bspline.c daubechies.c haar.c

TESTS = test

check_PROGRAMS = test

test_LDADD = libgslwavelet.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la  ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

test_SOURCES = test.c


