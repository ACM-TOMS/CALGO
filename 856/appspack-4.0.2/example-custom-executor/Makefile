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

OFILES = main.o master.o worker.o executor.o

appspack_custom: $(OFILES)
	$(CXX) -o appspack_custom $(OFILES) $(LDFLAGS) $(LIBS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o main.o main.cpp

master.o: master.cpp msgtags.hpp executor.hpp
	$(CXX) $(CXXFLAGS) -c -o master.o master.cpp

worker.o: worker.cpp msgtags.hpp
	$(CXX) $(CXXFLAGS) -c -o worker.o worker.cpp

executor.o: executor.cpp executor.hpp msgtags.hpp
	$(CXX) $(CXXFLAGS) -c -o executor.o executor.cpp

clean:
	-rm -f $(OFILES) appspack_custom
 
