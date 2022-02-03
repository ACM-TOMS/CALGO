# Makefile for the Level 1 BLAS test suite
# Requires a file etc/Makeinclude.arch for this architecture
# Some sample Makeinclude files are provided in the etc directory.
#
# The environment variable ARCH must be set to match, for example
#    setenv ARCH intel
# The object files, module files, executable files, and output
# files are all created in the ARCH directory by default.
#
# This organization allows you to test multiple architectures
# without overwriting the files created in previous tests.
#
# Edward Anderson
# May 6, 2016
#
ARCHFILE = etc/Makeinclude.${ARCH}
include $(ARCHFILE)
SRCDIR  = src
TSTDIR = testing

.SUFFIXES: .f .F .f90 .F90

all: $(ARCHFILE) lib test

# Exit if Makeinclude file does not exist
$(ARCHFILE):
	@echo "Check that the ARCH environment variable is set correctly"
	@exit 1

congen:
	( cd $(TSTDIR); make -f Makecongen check )

check: lib test
	( cd $(TSTDIR); make check )

lib:
	( cd $(SRCDIR); make )

test:
	( cd $(TSTDIR); make )

clean:
	( cd $(SRCDIR); make clean )
	( cd $(TSTDIR); make clean )
