# Define EPC F90 compiler and flags
#FC = epcf90
#FFLAGS = -C -d1 -g -temp=/tmp -u
#FFLAGS = -temp=/tmp -O

# Define Nag f90 compiler and flags
FC = f90
FFLAGS = -g

# Define rule for .f to .o and .f90 to .o
.SUFFIXES : .f .f90 .o
.f.o:
	$(FC) $(FFLAGS) -c $<
.f90.o:
	$(FC) $(FFLAGS) -c $<

all: res1 res2 res3 res4 res5 res6

res1: zmlib.o fmlib.o driver1.o
	$(FC) $(FFLAGS) zmlib.o fmlib.o driver1.o -o driver1
	driver1 > res1

res2: zmlib.o fmlib.o driver2.o
	$(FC) $(FFLAGS) zmlib.o fmlib.o driver2.o -o driver2
	driver2 > res2

res3: fmlib.o driver3.o
	$(FC) $(FFLAGS) fmlib.o driver3.o -o driver3
	driver3 > res3

res4: driver4.o fmlib.o
	$(FC) $(FFLAGS) driver4.o fmlib.o -o driver4
	driver4 > res4

res5: fmzmcomm.o fmzm90.o fmlib.o  zmlib.o driver5.o
	$(FC) $(FFLAGS) fmzmcomm.o fmzm90.o fmlib.o  zmlib.o driver5.o -o driver5
	driver5 > res5

res6: fmzmcomm.o fmzm90.o fmlib.o  zmlib.o driver6.o
	$(FC) $(FFLAGS) fmzmcomm.o fmzm90.o fmlib.o  zmlib.o driver6.o -o driver6
	driver6 > res6

clean: 
	rm  -rf driver4 driver6 driver3 driver5 driver1 driver2
	rm -rf  *.o *.LOG res*
