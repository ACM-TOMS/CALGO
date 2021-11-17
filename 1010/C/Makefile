export CC = gcc
export PERC = $(PWD)
export BIN = $(PERC)/bin
export SRC = $(PERC)/Src
export DRV = $(PERC)/Drivers
export OF = -O3
export CFLAGS=-Wall -Wextra

default: timingtest timingtest_sample_F accuracytest statanalysis simpletest 

COMOBJ=$(DRV)/Abramowitz.o $(DRV)/Flocke.o $(DRV)/NumRecipeHQRL.o $(DRV)/StrobachLDLT.o $(DRV)/common.o $(DRV)/FastQuarticSolver.o $(DRV)/NumRecipeHQR.o $(DRV)/Shmakov.o 

sources: 
	cd $(SRC) && $(MAKE)

solvers:
	cd $(DRV) && $(MAKE)

common: sources solvers 
	
accuracytest: common $(DRV)/accuracytest.c 
	$(CC) $(CFLAGS) $(OF) -o $(BIN)/accuracytest $(DRV)/accuracytest.c $(COMOBJ) $(SRC)/quartic_solver.o -lm

timingtest: common $(DRV)/timingtest.c 
	$(CC) $(CFLAGS) $(OF) -o $(BIN)/timingtest $(DRV)/timingtest.c $(COMOBJ) $(SRC)/quartic_solver.o -lm

timingtest_sample_F: common $(DRV)/timingtest_sample_F.c 
	$(CC) $(CFLAGS) $(OF) -o $(BIN)/timingtest_sample_F $(DRV)/timingtest_sample_F.c $(COMOBJ) $(SRC)/quartic_solver.o -lm

statanalysis: common $(DRV)/statanalysis.c
	$(CC) $(CFLAGS) $(OF) -o $(BIN)/statanalysis $(DRV)/statanalysis.c $(COMOBJ) $(SRC)/quartic_solver.o -lm
 
simpletest: sources $(DRV)/simpletest.c 
	$(CC) $(CFLAGS) -o $(BIN)/simpletest $(DRV)/simpletest.c $(SRC)/quartic_solver.o $(SRC)/quartic_solver_cmplx.o -lm

clean: 
	rm -f $(SRC)/*.o $(DRV)/*.o $(BIN)/simpletest $(BIN)/accuracytest $(BIN)/timingtest $(BIN)/statanalysis $(BIN)/timingtest_sample_F
        
