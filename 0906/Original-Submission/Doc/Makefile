#Generate library filename
LIB=libelrint3d.a
SRC=elrint3d.cc elrint3dSequence.cc precomputedCoords.cc
OBJ=$(SRC:.cc=.o)
AR= ar cq
CXXFLAGS=-Wall -Weffc++ -pedantic
CXXFLAGS=-Wall -Weffc++ -O3
CPPFLAGS=-I. 
LDFLAGS=-L. -lelrint3d
COMPILE.cc=$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(TARGET_ARCH)

#Make binary executable for Examples 1-5
all:    $(LIB) 
	$(COMPILE.cc) -o test01 test01.cc   $(LDFLAGS)
	$(COMPILE.cc) -o test02 test02.cc   $(LDFLAGS)
	$(COMPILE.cc) -o test03 test03.cc   $(LDFLAGS)
	$(COMPILE.cc) -o test04 test04.cc   $(LDFLAGS)
	$(COMPILE.cc) -o test05 test05.cc   $(LDFLAGS)	

precomputedCoords.o: precomputedCoords.cc precomputedCoords.h
	$(COMPILE.cc) -c $<

elrint3d.o: elrint3d.cc elrint3d.h cubUtil.h infiniteMap.h elrint3dSequence.h\
            embeddedSequence.h embeddedCubRule.h integrand.h
	$(COMPILE.cc) -c  $<
	
elrint3dSequence.o: elrint3dSequence.cc elrint3dSequence.h precomputedCoords.h\
                    embeddedSequence.h		    
	$(COMPILE.cc) -c $<

#Make a library 
$(LIB): $(OBJ)
	$(AR) $@ $(OBJ)
	ranlib $@
	
lib: $(LIB)

#Clean the package	
clean: 
	rm -rf *.o $(LIB)  test01 test02 test03 test04 test05
