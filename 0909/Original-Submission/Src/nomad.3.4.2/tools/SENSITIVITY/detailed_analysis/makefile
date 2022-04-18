EXE                 = detailed_analysis.exe
COMPILATOR          = g++
COMPILATOR_OPTIONS  = -ansi -Wall -O3 -DGCC_X
L1                  = $(NOMAD_HOME)/lib/nomad.a
LIBS                = $(L1) -lc -lm
INCLUDE             = -I$(NOMAD_HOME)/src -I.
COMPILE             = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
OBJS                = detailed_analysis.o

$(EXE): $(OBJS)
	$(COMPILATOR) -o $(EXE) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS)

detailed_analysis.o: detailed_analysis.cpp $(L1)
	$(COMPILE) detailed_analysis.cpp

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
