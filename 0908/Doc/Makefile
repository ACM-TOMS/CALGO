%.o : %.cpp
	$(CPP) $(CPPOPTS) -c $*.cpp

# For SPARC Solaris with Sun C++ compiler
CPP = CC
CPPOPTS = -O1 -DREV

# For SPARC Solaris and PowerPC with GNU C++ compiler
#CPP = g++
#CPPOPTS = -O1 -DREV

# For x86 Linux and Cygwin Windows with GNU C++ compiler
#CPP = g++
#CPPOPTS = -O1 -DDOUBLE

# For x86 with Mac OS X
#CPP = g++
#CPPOPTS = -O1


OBJS = main.o ExactSum.o

res: driver
	driver 1000 10 1 > res
	driver 1000 30 2 >> res
	driver 1000 50 3 >> res
	driver 1000 1000 4 >> res

driver: $(OBJS)
	$(CPP) $(CPPOPTS) -o driver $(OBJS)

clean:
	rm -rf main *.o

main.o: main.cpp ExactSum.h
	$(CPP) $(CPPOPTS) -c main.cpp

ExactSum.o: ExactSum.cpp ExactSum.h
	$(CPP) $(CPPOPTS) -c ExactSum.cpp
