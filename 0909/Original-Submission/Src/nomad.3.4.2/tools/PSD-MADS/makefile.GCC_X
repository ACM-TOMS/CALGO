EXE                 = psdmads
COMPILATOR          = mpic++
COMPILATOR_OPTIONS  = -ansi -Wall -O3
L1                  = $(NOMAD_HOME)/lib/nomad.a
LIBS                = $(L1) -lc -lm
INCLUDE             = -I$(NOMAD_HOME)/src -I.
COMPILE             = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
OBJS                = main.o Cache_Server.o Master_Slaves.o

$(EXE): $(OBJS)
	$(COMPILATOR) -o $(EXE) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS)

main.o: main.cpp Master_Slaves.o
	$(COMPILE) main.cpp

Master_Slaves.o: Master_Slaves.cpp Master_Slaves.hpp Cache_Server.o
	$(COMPILE) Master_Slaves.cpp

Cache_Server.o: Cache_Server.cpp Cache_Server.hpp $(L1)
	$(COMPILE) Cache_Server.cpp

clean:
	@echo "cleaning obj files"
	@rm -f $(OBJS)

del:
	@echo "cleaning trash files"
	@rm -f core *~
	@echo "cleaning obj files"
	@rm -f $(OBJS)
	@echo "cleaning exe file"
	@rm -f $(EXE)