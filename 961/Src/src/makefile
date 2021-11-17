####################################################################
#  SHHEIG routines makefile                                        #
#  Makefile for creating/updating the SHHEIG Library object file   # 
#  on Unix machines.                                               #
#  SHHEIG, Release 1.0                      ./shheig/src/makefile  #
#  Matthias Voigt, MPI Magdeburg                                   #
#  July 23, 2013.                                                  #
####################################################################
#
#  This is the makefile to create/update the library for SHHEIG.
#  The SHHEIG Library routines are written for double precision only.
#
#  The command
#       make
#  without any arguments creates or updates a library called
#       slicot.a
#  in the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make double FRC=FRC
#
#######################################################################

include ../make.inc

DSLSRC = \
    DGHFDF.o DGHFET.o DGHFEX.o DGHFEY.o DGHFST.o DGHFXC.o DGHFYR.o \
    DGHUDF.o DGHUDP.o DGHUET.o DGHUEX.o DGHUEY.o DGHURV.o DGHUSP.o \
    DGHUST.o DGHUTP.o DGHUTR.o DGHUXC.o DGHUXP.o DGHUYR.o DLACPV.o \
    ZGHFDF.o ZGHFEX.o ZGHFEY.o ZGHFXC.o ZGHUDF.o ZGHUDP.o ZGHUEX.o \
    ZGHUEY.o ZGHUXC.o ZGHUXP.o ZLACPV.o

all: double

double: $(DSLSRC)
	$(ARCH) $(ARCHFLAGS) $(SHHEIGLIB) $(DSLSRC)

$(DSLSRC): $(FRC)

FRC:
	@FRC=$(FRC)

clean:
	rm -f *.o

.f.o: 
	$(FORTRAN) $(OPTS) -c $<
