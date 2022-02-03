EXE                 = nomad_for_dll.exe
COMPILATOR          = g++
COMPILATOR_OPTIONS  = -O3 -DGCC_WINDOWS
L1                  = $(NOMAD_HOME)/lib/nomad.lib
LIBS                = $(L1) -lm
INCLUDE             = -I$(NOMAD_HOME)/src -I.
COMPILE             = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
OBJS                = nomad_for_dll.o

$(EXE): $(OBJS)
	$(COMPILATOR) -o $(EXE) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS)

mpi: EXE                 = nomad_for_dll_MPI.exe
mpi: COMPILATOR_OPTIONS += -DUSE_MPI
mpi: INCLUDE            += -I$(MPI_HOME)/include
mpi: L1                  = $(NOMAD_HOME)/lib/nomad_MPI.lib
mpi: LIBS                = $(L1) -lm $(MPI_HOME)/lib/mpi.lib
mpi: clean $(OBJS)
	@echo "   building the MPI version ..."
	@echo "   exe file : "$(EXE)
	@$(COMPILATOR) -o  $(EXE) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS)
	@strip  $(EXE)

nomad_for_dll.o: nomad_for_dll.cpp $(L1)
	$(COMPILE) nomad_for_dll.cpp

clean:
	@echo cleaning obj files
	@del $(OBJS)

del:
	@echo cleaning obj files
	@del $(OBJS)
	@echo cleaning exe file
	@del $(EXE)
