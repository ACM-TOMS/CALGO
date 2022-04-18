# Set the location of the APPSPACK hearder files
APPSPACKINCDIR = ../src

# Set the location of the APPSPACK libraries
APPSPACKLIBDIR = ../src

# Set the C++ compiler. 
# This must be the same as what was used to compile APPSPACK. 
# To check which compiler was used, in the APPSPACK directory, 
# execute: grep ^"CXX =" Makefile
CXX = mpiCC

CXXFLAGS = -I$(APPSPACKINCDIR) -DHAVE_CONFIG_H
LDFLAGS = -L$(APPSPACKLIBDIR)
LIBS = -lappspack

appspack_custom: main.o evaluator.o
	$(CXX) -o appspack_custom main.o evaluator.o $(LDFLAGS) $(LIBS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o main.o main.cpp

evaluator.o: evaluator.cpp evaluator.hpp
	$(CXX) $(CXXFLAGS) -c -o evaluator.o evaluator.cpp

clean:
	-rm main.o evaluator.o appspack_custom
