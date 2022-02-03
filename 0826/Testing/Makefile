############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          Makefile
#
#  Purpose:         Eigenroutine Testing Makefile
#
#  Creation date:   March 20, 1995
#
#  Modified:
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../../SLmake.inc

F77 = f90
F77LOADER = f90

znepexe = xznep
zevcexe = xzevc

zmatgen = pzmatgen.o pmatgeninc.o

znep =  pznepdriver.o pznepinfo.o pznepfchk.o $(zmatgen)
zevc =  pzevcdriver.o pzevcinfo.o pzget22.o $(zmatgen) 

all : complex16

complex: $(cnepexe) $(cevcexe)
complex16: $(znepexe) $(zevcexe)

$(cnepexe) : $(cnep)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cnepexe) $(cnep) $(PARLIBS)

$(cevcexe) : $(cevc)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cevcexe) $(cevc) $(PARLIBS)

$(znepexe) : $(znep)
	$(F77LOADER) $(F77LOADFLAGS) -o $(znepexe) $(znep) $(PARLIBS)

$(zevcexe) : $(zevc)
	$(F77LOADER) -g $(F77LOADFLAGS) -o $(zevcexe) $(zevc) $(PARLIBS)

$(cnep): $(FRC)
$(cevc): $(FRC)
$(znep): $(FRC)
$(zevc): $(FRC)

FRC:
	@FRC=$(FRC)

clean :
	rm -f *.o

.c.o:
	$(CC) $(CDEFS) $(OPTS) -c $*.c

.f.o: 
	$(F77) -c $(OPTS) $*.f
