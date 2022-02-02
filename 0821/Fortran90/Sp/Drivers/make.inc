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
all: test1 test2 test3 test4 bench1

test1: $(FLIBS) test1.o
	$(FC) $(FLOADFLAGS) -o test1 test1.o $(LIBS)

test2: $(FLIBS) test2.o
	$(FC) $(FLOADFLAGS) -o test2 test2.o $(LIBS)

test3: $(FLIBS) test3.o
	$(FC) $(FLOADFLAGS) -o test3 test3.o $(LIBS)

test4: $(FLIBS) test4.o
	$(FC) $(FLOADFLAGS) -o test4 test4.o $(LIBS)

bench1: $(FLIBS) bench1.o
	$(FC) $(FLOADFLAGS) -o bench1 bench1.o $(LIBS)

clean:
	rm -f *.o test1 test2 test3 test4 bench1 *.$(MODULES)
