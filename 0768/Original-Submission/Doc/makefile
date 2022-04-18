#
#   TENSOLVE directory
#

FFLAGS = -O -u

#   Files for TENSOLVE
#   The program runs nonlinear equations and nonlinear
#   least squares problems

EXAMPLE1 = driver1.o
EXAMPLE2 = driver2.o
TENSOLVE = tensolve.o
UNCMIN = uncmin.o
BLAS = blas.o

ex1 : $(EXAMPLE1) $(TENSOLVE) $(UNCMIN) $(BLAS)
	f77 $(FFLAGS) $(EXAMPLE1) $(TENSOLVE) $(UNCMIN) \
	$(BLAS) -o tensolve
ex2 : $(EXAMPLE2) $(TENSOLVE) $(UNCMIN) $(BLAS)
	f77 $(FFLAGS) $(EXAMPLE2) $(TENSOLVE) $(UNCMIN) \
	$(BLAS) -o tensolve

