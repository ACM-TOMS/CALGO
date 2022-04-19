# Uses a simple non-recursive make, where all makefiles are in or directly below
# the top directory
#
# Making from the blas-rmd top directory:
#   'make' makes the rmd-library lib/libblasrmd.a and makes and runs demos
#   'make src' just makes the rmd-library
#   'make all' makes everything (including the reference BLAS) and runs tests and demos
#   'make allbutdoc' make all but the documentation
#   'make testrun' makes and runs all tests. No output means that test passed.
#   'make timerun' makes and runs the timing program with default settings
#   'make doc' remakes the documentation
#   'make clean' removes everything makeable except pdf-documentation
#   'make clean-all' removes also pdf-doc
#   'make <directory>' makes the default target of <directory> and all its dependencies
#   'make <target>' makes a specified target (in a directory) and all its dependencies
#
# Making in a subdirectory:
#   'make' makes the default target of the directory, assuming its dependencies are made
#   'make <target>' makes the specified target
#
# Make options:
#   Two make parameters can be used to set compile switches:
#     make debug=1 [target]     Use debug switches
#     make optimize=1 [target]  Maximum optimization

#   The compiler and additional options may be set with:
#     make FC=<fortrancompiler> [target]
#     make FOPTS=<options> [target]

#   The default is to make with all options off (-O0). For timing, use:
#     make clean; make time optimize=1
#
# Environment variables:
#     $BLAS can be set to openblas, mkl, accelerate or refblas, e.g. in .bashrc
#     $OPENBLAS can be set to the installation root of openblas, in .bashrc
#     $FC can be set to the Fortran compiler
#     $FOPTS can be set to additional compile options

#   These variables can also be changed in 'make.inc' or set on the make command
#   line, e.g. 'make demo BLAS=mkl'

TOPLEVEL:=true
include make.inc

ALLDIR := src tools demo test time doc refblas article
default: demorun
all: $(ALLDIR) demorun testdemo testrun
allbutdoc: $(filter-out doc article,$(ALLDIR))

$(foreach DIR,$(ALLDIR),$(eval include $(DIR)/Makefile))

clean::
	rm -f make.log *.mod
