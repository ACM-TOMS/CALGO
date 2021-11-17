# Makefile to create a main program that tests RELIADIFF on a laplace Transforms Database
# Version 1.0


ifneq ($(CC),g++)
CC = g++
endif

ifndef rootdir 
rootdir = ../../
endif

ifndef RELdir
RELdir = $(rootdir)/src
endif

ifndef incdir 
incdir = $(rootdir)
endif


ifndef dbdir
dbdir = ./database
endif

ifndef utildir
utildir = ../utility
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

ifndef gsldirlib
gsldirlib = $(gsllib)/lib
endif

ifndef CFLAGS
CFLAGS = -Wall -ansi -pedantic -g
endif

ifndef RM
RM = rm -f -r
endif


LIBS = -lm -L$(RELdir) -lreliadiff -L$(dbdir) -ldatabase -L$(utildir) -lutil  -L$(gsldirlib) -lgsl -lgslcblas 

OBJ = TEST_ON_DATABASE.o 
BIN = test_on_database

all: $(BIN) 

relib: 
	@echo
	@echo
	@echo "############################################################################"
	@echo "RELIADIFF DEMO: creating RELIADIFF library"
	@echo "############################################################################"
	@$(MAKE) -C $(RELdir)

util: $(gsldirinc)/gsl relib
	@echo
	@echo
	@echo "############################################################################"
	@echo "RELIADIFF DEMO: creating utility library"
	@echo "############################################################################"
	@$(MAKE) -C $(utildir) 
	
db: util 
	@echo
	@echo
	@echo "############################################################################"
	@echo "RELIADIFF DEMO: creating Transform/Inverse Databases library"
	@echo "############################################################################"
	@$(MAKE) -C $(dbdir) 

$(OBJ): TEST_ON_DATABASE.c db
	@echo
	@echo
	@echo "############################################################################"
	@echo "RELIADIFF DEMO: compiling main program to test RELIADIFF on database..."
	@echo "############################################################################"
	@$(CC) $(CFLAGS) -I$(gsldirinc) -I$(incdir) -c TEST_ON_DATABASE.c 
	
$(BIN): $(OBJ) 
	@echo
	@echo
	@echo "############################################################################"
	@echo "RELIADIFF DEMO: linking main program to test RELIADIFF on database..."
	@echo "############################################################################"
	@$(CC) $(CFLAGS) -o $(BIN) $(OBJ) $(LIBS)
	@${RM} $(OBJ) 
	@echo
	@echo
	@echo "############################################################################"
	@echo "############################################################################"
	@echo "############################################################################"
	@echo "############################################################################"
	@echo "Test driver installed!!!"
	@echo "Now you can run the program $(BIN) or a script calling it in this directory."
	@echo "(Consult the RELIADIFF_userguide.pdf file to know about the arguments)"
	@echo "############################################################################"
	@echo "############################################################################"
	@echo "############################################################################"
	@echo "############################################################################"

clean: 
	@echo "RELIADIFF TEST_ON_DATABASE: deleting existing executable and output files..."
	@${RM} $(BIN) *.txt $(dbdir)/*.a $(RELdir)/*.a $(utildir)/*.a







	
	
