EXE                 = test.exe
COMPILATOR          = g++
FCOMPILATOR         = gfortran
COMPILATOR_OPTIONS  = -ansi -Wall -O3 -DGCC_X
L1                  = $(NOMAD_HOME)/lib/nomad.a
LIBS                = $(L1) -lc -lm
INCLUDE             = -I$(NOMAD_HOME)/src -I.
COMPILE             = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
FCOMPILE            = $(FCOMPILATOR) -O3 -c
OBJS                = nomad.o test.o

$(EXE): $(OBJS)
	$(FCOMPILATOR) -o $(EXE) $(OBJS) $(LIBS) -lstdc++

nomad.o: nomad.cpp $(L1)
	$(COMPILE) nomad.cpp

test.o: test.f
	$(FCOMPILE) test.f

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
