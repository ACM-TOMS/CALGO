####################################################################
#  SHHEIG examples makefile                                        #
#  Makefile for generating and running SHHEIG Library example      # 
#  programs on Unix machines.                                      #
#  SHHEIG, Release 1.0                 ./slicot/examples/makefile  #
#  Matthias Voigt, MPI Magdeburg                                   #
#  July 23, 2013.                                                  #
####################################################################
#
#  This makefile compiles, links, and runs the example programs for 
#  the SHHEIG Library on Unix machines.
#
#  The example programs can be executed for double precision only.
#  To compile, link, and run the example programs, enter 'make'.
#  The executable files are created in the current directory level.
#  The files with the results have the extension .exa, and are also
#  created in the current directory level, so they can automatically
#  be compared with the .res files provided in this directory. Note
#  that, for some programs, the signs of some matrix elements could
#  differ; this does not mean erroneous results.
#
#  To remove the .exa files as well as the object files and
#  excutable programs, enter
#       make clean
#
#  To re-run specific programs after a make, enter (for example):
#       'rm TDGHFDF; make'  or:
#       'make TDGHFDF' or:
#       'touch DGHFDF.dat; make' (to re-run the TDGHFDF program).
#
####################################################################

include ../make.inc

all: DGHFDF.exa DGHFST.exa DGHUDF.exa DGHURV.exa DGHUST.exa DGHUTR.exa \
     ZGHFDF.exa ZGHUDF.exa

DGHFDF.exa: DGHFDF.dat TDGHFDF; ./TDGHFDF <DGHFDF.dat >DGHFDF.exa
DGHFST.exa: DGHFST.dat TDGHFST; ./TDGHFST <DGHFST.dat >DGHFST.exa
DGHUDF.exa: DGHUDF.dat TDGHUDF; ./TDGHUDF <DGHUDF.dat >DGHUDF.exa
DGHURV.exa: DGHURV.dat TDGHURV; ./TDGHURV <DGHURV.dat >DGHURV.exa
DGHUST.exa: DGHUST.dat TDGHUST; ./TDGHUST <DGHUST.dat >DGHUST.exa
DGHUTR.exa: DGHUTR.dat TDGHUTR; ./TDGHUTR <DGHUTR.dat >DGHUTR.exa
ZGHFDF.exa: ZGHFDF.dat TZGHFDF; ./TZGHFDF <ZGHFDF.dat >ZGHFDF.exa
ZGHUDF.exa: ZGHUDF.dat TZGHUDF; ./TZGHUDF <ZGHUDF.dat >ZGHUDF.exa

TDGHFDF: TDGHFDF.o; $(LOADER) -o TDGHFDF TDGHFDF.o $(LOADOPTS)
TDGHFST: TDGHFST.o; $(LOADER) -o TDGHFST TDGHFST.o $(LOADOPTS)
TDGHUDF: TDGHUDF.o; $(LOADER) -o TDGHUDF TDGHUDF.o $(LOADOPTS)
TDGHURV: TDGHURV.o; $(LOADER) -o TDGHURV TDGHURV.o $(LOADOPTS)
TDGHUST: TDGHUST.o; $(LOADER) -o TDGHUST TDGHUST.o $(LOADOPTS)
TDGHUTR: TDGHUTR.o; $(LOADER) -o TDGHUTR TDGHUTR.o $(LOADOPTS)
TZGHFDF: TZGHFDF.o; $(LOADER) -o TZGHFDF TZGHFDF.o $(LOADOPTS)
TZGHUDF: TZGHUDF.o; $(LOADER) -o TZGHUDF TZGHUDF.o $(LOADOPTS)

$(TDGHFDF): $(FRC)
$(TDGHFST): $(FRC)
$(TDGHUDF): $(FRC)
$(TDGHURV): $(FRC)
$(TDGHUST): $(FRC)
$(TDGHUTR): $(FRC)
$(TZGHFDF): $(FRC)
$(TZGHUDF): $(FRC)

FRC:
   @FRC=$(FRC)

clean:
	rm -f *.exa *.o \
     TDGHFDF TDGHFST TDGHUDF TDGHURV TDGHUST TDGHUTR TZGHFDF TZGHUDF

TDGHFDF.o: TDGHFDF.f; $(FORTRAN) $(OPTS) -c $<
TDGHFST.o: TDGHFST.f; $(FORTRAN) $(OPTS) -c $<
TDGHUDF.o: TDGHUDF.f; $(FORTRAN) $(OPTS) -c $<
TDGHURV.o: TDGHURV.f; $(FORTRAN) $(OPTS) -c $<
TDGHUST.o: TDGHUST.f; $(FORTRAN) $(OPTS) -c $<
TDGHUTR.o: TDGHUTR.f; $(FORTRAN) $(OPTS) -c $<
TZGHFDF.o: TZGHFDF.f; $(FORTRAN) $(OPTS) -c $<
TZGHUDF.o: TZGHUDF.f; $(FORTRAN) $(OPTS) -c $<

.f.o: ; $(FORTRAN) $(OPTS) -c $<

