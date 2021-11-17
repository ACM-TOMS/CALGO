# Makefile to create a Laplace Transforms Database
# Version 1.0

ifneq ($(CC),g++)
CC = g++
endif

ifndef rootdir 
rootdir = ../../..
endif

ifndef utildir
utildir = ../../utility
endif

ifndef incdir 
incdir = $(rootdir)
endif

ifndef dirlib
dirlib  = .
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
RM = rm -f
endif


OBJ = dbInvLaplace.o dbLaplace.o
BIN = libdatabase.a

all: $(BIN) 
	
dbInvLaplace.o: dbInvLaplace.c dbL.h $(utildir)/Util.h $(incdir)/src/RELIADIFF.h
	@$(CC) $(CFLAGS) -I$(gsldirinc) -I$(incdir) -c dbInvLaplace.c 
dbLaplace.o: dbLaplace.c dbL.h $(utildir)/Util.h $(incdir)/src/RELIADIFF.h
	@$(CC) $(CFLAGS) -I$(gsldirinc) -I$(incdir) -c dbLaplace.c 

$(BIN): $(OBJ)
	@echo "RELIADIFF: compiling functions database..."
	@ar qc $(BIN) $(OBJ)
	@${RM} $(OBJ)
	
clean: 
	@echo "RELIADIFF: deleting the built database file..."
	@${RM} $(BIN)


	




