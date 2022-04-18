EXE                 = test.exe
COMPILATOR          = g++
FCOMPILATOR         = gfortran
COMPILATOR_OPTIONS  = -ansi -Wall -O3 -DGCC_X
L1                  = $(NOMAD_HOME)/lib/nomad.a
LIBS                = $(L1) -lc -lm
INCLUDE             = -I$(NOMAD_HOME)/src -I.
COMPILE             = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
FCOMPILE            = $(FCOMPILATOR) -O3 -c
OBJS                = test.o bb.o

$(EXE): $(OBJS)
	$(COMPILATOR) -o $(EXE) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS) -lgfortran

test.o: test.cpp bb.o $(L1)
	$(COMPILE) test.cpp

bb.o: bb.f
	$(FCOMPILE) bb.f

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
