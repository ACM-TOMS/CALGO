# Makefile of RELIADIFF Utility Functions
# Version 1.0

ifneq ($(CC),g++)
CC = g++
endif

ifndef rootdir 
rootdir = ../../
endif

ifndef incdir 
incdir = $(rootdir)
endif

ifndef dirlib
dirlib = .
endif

## CHECK
#############################################################################################################
### GSL building directory: change it if you know a different one in your envirmoment                     ###	
ifndef gsllib																							  
gsllib = /usr/local
endif
#############################################################################################################

ifndef gsldirinc
gsldirinc = $(gsllib)/include
endif


ifndef CFLAGS
CFLAGS = -Wall -ansi -pedantic -g
endif

ifndef RM
RM = rm -f
endif

OBJU = Util.o 
BINU = libutil.a



all: $(BINU)
	
Util.o: Util.c Util.h $(gsldirinc)/gsl/ 
	@echo "RELIADIFF: compiling utility functions..."
	@$(CC) $(CFLAGS) -I$(gsldirinc) -I$(incdir) -c Util.c

$(BINU): $(OBJU)
	@echo "RELIADIFF: creating utility library..."
	@ar qc $(BINU) $(OBJU)
	@${RM} $(OBJU) 
	
clean: 
	@echo "RELIADIFF: deleting built utility files..."
	@${RM} $(dirlib)/$(BINU)

