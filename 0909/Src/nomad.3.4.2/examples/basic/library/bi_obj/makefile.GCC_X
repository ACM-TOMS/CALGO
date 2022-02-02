EXE                 = basic_lib
COMPILATOR          = g++
COMPILATOR_OPTIONS  = -ansi -Wall -O3 -DGCC_X
L1                  = $(NOMAD_HOME)/lib/nomad.a
LIBS                = $(L1) -lc -lm
INCLUDE             = -I$(NOMAD_HOME)/src -I.
COMPILE             = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
OBJS                = basic_lib.o

$(EXE): $(OBJS)
	$(COMPILATOR) -o $(EXE) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS)

basic_lib.o: basic_lib.cpp $(L1)
	$(COMPILE) basic_lib.cpp

mpi: COMPILATOR          = mpic++
mpi: EXE                 = basic_lib_MPI
mpi: L1                  = $(NOMAD_HOME)/lib/nomad.MPI.a
mpi: LIBS                = $(L1) -lc -lm
mpi: COMPILATOR_OPTIONS += -DUSE_MPI
mpi: clean $(OBJS)
	@echo "   building the MPI version ..."
	@echo "   exe file : "$(EXE)
	@$(COMPILATOR) -o  $(EXE) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS)
	@strip  $(EXE)

clean:
	@echo "   cleaning obj files"
	@rm -f $(OBJS)

del:
	@echo "   cleaning trash files"
	@rm -f core *~
	@echo "   cleaning obj files"
	@rm -f $(OBJS)
	@echo "   cleaning exe file"
	@rm -f $(EXE)
