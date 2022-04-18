# Makefile of RELIADIFF
# Version 1.1

# This Makefile creates a library to link to when compiling a code using RELIADIFF function

ifneq ($(CC),g++)
CC = g++
endif

ifndef rootdir 
rootdir = ..
endif

ifndef incdir
incdir = $(rootdir)
endif

ifndef RELdir 
RELdir = .
endif

ifndef dirlib
dirlib = .
endif

ifndef CFLAGS
CFLAGS = -Wall -ansi -pedantic -g
endif

ifndef RM
RM = rm -f
endif

OBJ = RELIADIFF.o phi.o
BIN = libreliadiff.a


all: $(BIN)
	
phi.o: phi.c RELIADIFF.h fadbad/tadiff.h
	@echo "RELIADIFF: compiling phi..."
	@$(CC) $(CFLAGS) -I$(incdir) -c phi.c 
		
RELIADIFF.o: RELIADIFF.c phi.c RELIADIFF.h fadbad/tadiff.h
	@echo "RELIADIFF: compiling RELIADIFF..."
	@$(CC) $(CFLAGS) -I$(incdir) -c RELIADIFF.c 

$(BIN): $(OBJ)
	@echo "RELIADIFF: creating library..."
	@ar qc $(BIN) $(OBJ)
	@${RM} $(OBJ) 
	
clean: 
	@echo "RELIADIFF: deleting the built library..."
	@${RM} $(dirlib)/$(BIN) *.o



