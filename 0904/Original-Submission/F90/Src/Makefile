##############################################################################
############################################################################## 
##                                                                          ##
##             Makefile for building the SCASY library                      ##
##                                                                          ##
##############################################################################
##############################################################################
#
include ./../Make_include
#
all:		libscasy
#
OBJECTS = BDLAAPP.o DLAQR1.o DTRGSYL.o DTRSYUP.o PDEXTCHK.o PDHSEQR.o \
          PDLAQR1.o PDLAQR5.o PGEGCSYD.o PGESYCTD.o PIPARMQ.o PSYCTCON.o \
          PTRLYCTD.o BDLAEXC.o DLAQR6.o DTRLYCT.o DUBEXMA.o PDGEBAL.o \
          PDIMPRED.o PDMATGEN2.o PGEGLYCTD.o PGESYDTD.o PLYCTCON.o \
          PSYDTCON.o PTRLYDTD.o BDTREXC.o DLARAN.o DTRLYDT.o P1SQMATGD.o \
          PDGEES.o PDLAHQR.o PDLAQR3.o PDORMHR.o PGEGLYDTD.o PGLYCTCON.o \
          PLYDTCON.o PTRGCSYD.o PTRSYCTD.o DBEXMAT.o DLAREF.o DTRMX.o \
          P2SQMATGD.o PDGEHRD.o PDLAMVE.o PDQZINC.o PGEGSYLD.o \
          PGLYDTCON.o PMATGENINC.o PTRGLYCTD.o PTRSYDTD.o DGEMX.o DSECND.o \
          DTRSYCT.o PBDTRSEN.o PDGGHRD.o PDLAPRNT.o PDROT.o \
          PGELYCTD.o PGSYLCON.o PQRHLPS.o PTRGLYDTD.o SCASYHLPS.o DLAPRNT.o \
          DTRGCSY.o DTRSYDT.o PDBCKRD.o PDHGEQZ.o PDLAQR0.o \
          PGCSYCON.o PGELYDTD.o PILAENVX.o PQRRMMM.o PTRGSYLD.o PBDTGSEN.o \
          PSB02MD.o PSB02MU.o \
          limits.o core-unlimited.o core-unlimited_.o
#
FILES   = BDLAAPP.f DLAQR1.f DTRGSYL.f DTRSYUP.f PDEXTCHK.f PDHSEQR.f \
          PDLAQR1.f PDLAQR5.f PGEGCSYD.f PGESYCTD.f PIPARMQ.f PSYCTCON.f \
          PTRLYCTD.f BDLAEXC.f DLAQR6.f DTRLYCT.f DUBEXMA.f PDGEBAL.f \
          PDIMPRED.f PDMATGEN2.f PGEGLYCTD.f PGESYDTD.f PLYCTCON.f \
          PSYDTCON.f PTRLYDTD.f BDTREXC.f DLARAN.f DTRLYDT.f P1SQMATGD.f \
          PDGEES.f PDLAHQR.f PDLAQR3.f PDORMHR.f PGEGLYDTD.f PGLYCTCON.f \
          PLYDTCON.f PTRGCSYD.f PTRSYCTD.f DBEXMAT.f DLAREF.f DTRMX.f \
          P2SQMATGD.f PDGEHRD.f PDLAMVE.f PDQZINC.f PGEGSYLD.f \
          PGLYDTCON.f PMATGENINC.f PTRGLYCTD.f PTRSYDTD.f DGEMX.f DSECND.f \
          DTRSYCT.f PBDTRSEN.f PDGGHRD.f PDLAPRNT.f PDROT.f \
          PGELYCTD.f PGSYLCON.f PQRHLPS.f PTRGLYDTD.f SCASYHLPS.f DLAPRNT.f \
          DTRGCSY.f DTRSYDT.f PDBCKRD.f PDHGEQZ.f PDLAQR0.f \
          PGCSYCON.f PGELYDTD.f PILAENVX.f PQRRMMM.f PTRGSYLD.f PBDTGSEN.f \
	  PSB02MD.f PSB02MU.f

#
libscasy:	${OBJECTS}
		ar r libscasy.a ${OBJECTS} 
#
$(OBJECTS):	$(FILES) limits.c core-unlimited.c
		$(FC) $(FFLAGS) $(INCLUDES) -c $(FILES)
		$(CC) -c limits.c core-unlimited.c core-unlimited_.c
#
clean:		
		make cleanlib
		make cleanobj
#
cleanlib:
		$(RM) $(TARGET)
#
cleanobj:	
		$(RM) $(OBJECTS)
#
##############################################################################
##############################################################################

