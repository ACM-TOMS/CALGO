ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = foreign
if HAVE_MKL
  MKL_EXAMPLES = examples/implicit_mkl examples/brusselator_mkl examples/brusselator_radau_mkl
endif
if HAVE_GSL
  GSL_EXAMPLES = examples/brusselator_gsl
endif


SUBDIRS = src . examples/explicit examples/implicit $(MKL_EXAMPLES) $(GSL_EXAMPLES)

# DOXYGEN SUPPORT
include $(top_srcdir)/aminclude.am

# ensure the distribution of the doxygen configuration file
EXTRA_DIST = Doxyfile.in \
	examples/explicit/ref \
	examples/implicit/ref \
	examples/implicit_mkl/ref  \
	examples/brusselator_mkl/ref \
	examples/brusselator_gsl/ref \
	examples/explicit/convergence.m \
	examples/implicit/convergence.m \
	examples/implicit_mkl/convergence.m \
	examples/brusselator_mkl/convergence.m \
	examples/brusselator_gsl/convergence.m \
	examples/explicit/run.sh \
	examples/implicit/run.sh \
	examples/implicit_mkl/run.sh \
	examples/brusselator_mkl/run.sh \
	examples/brusselator_gsl/run.sh \
	examples/brusselator_mkl/timing.sh \
	examples/brusselator_radau_mkl/timing.sh \
	examples/implicit_mkl/implicit.h \
	examples/brusselator_mkl/brusselator.h\
	examples/brusselator_gsl/brusselator.h \
	doxygen-doc
