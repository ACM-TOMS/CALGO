#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!                                                                          !!!
#!!  This program is furnished by the U.S. Army Engineer Research and        !!!
#!!  Development Center, Major Shared Resource Center (ERDC MSRC) "as is"    !!!
#!!  and is accepted and used by the recipient with the express              !!!
#!!  understanding that the Government makes no warranties, expressed or     !!!
#!!  implied, concerning the accuracy, completeness, reliability, usability  !!!
#!!  or suitability for any particular purpose of the information and data   !!!
#!!  within this program or furnished in connection therewith, and the       !!!
#!!  Government shall be under no liability whatsoever to any person by      !!!
#!!  reason of any use made thereof. This program belongs to the U.S.        !!!
#!!  Government; therefore, the recipient further agrees not to assert any   !!!
#!!  proprietary rights therein or to represent the source code to anyone    !!!
#!!  as other than a Government program.                                     !!!
#!!                                                                          !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!  Authors:                                                                !!!
#!!                                                                          !!!
#!!  Richard J. Hanson (koolhans@rice.edu)                                   !!!
#!!  Rice University,  Rice Center for High Performance Software Research    !!!
#!!                                                                          !!!
#!!  Clay P. Breshears (clay.breshears@intel.com)                            !!!
#!!  KAI Software, a division of Intel Americas, Inc.                        !!!
#!!                                                                          !!!
#!!  Henry A. Gabb (henry.gabb@intel.com)                                    !!!
#!!  KAI Software, a division of Intel Americas, Inc.                        !!!
#!!                                                                          !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

OBJS = fpthrd.o ptf90.o

all: $(OBJS)

fpthrd.o:  fpthrd.f
	rm -f fpthrd.f90
#	cpp $(FPP) fpthrd.f fpthrd.f90
	cp fpthrd.f fpthrd.f90
	$(FC) $(FFLAGS) -c fpthrd.f90
	mv *.o *.$(MODULES) $(FLIBPATH)

ptf90.o:  ptf90.c
	$(CC) $(CFLAGS) -c ptf90.c
	mv *.o $(FLIBPATH)

clean:
	(cd $(FLIBPATH); rm -f $(OBJS) *.$(MODULES) core)
