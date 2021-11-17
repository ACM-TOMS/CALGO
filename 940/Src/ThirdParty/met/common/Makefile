# -*- Makefile -*-

include ${HOME}/makehead/make$(EXT)

all : mat2$(EXT).o mat3$(EXT).o testdoublemin$(EXT).out testrand$(EXT).out

mat2$(EXT).o: mat2.cpp mat2.h
	$(CPP) -o mat2$(EXT).o -c mat2.cpp

mat3$(EXT).o: mat3.cpp mat3.h
	$(CPP) -o mat3$(EXT).o -c mat3.cpp

testdoublemin$(EXT) : testdoublemin.C
	$(CPP) -o $@ testdoublemin.C $(LINKOPT)

testdoublemin$(EXT).out : testdoublemin$(EXT)
	$? > $@

testrand$(EXT) : testrand.C
	$(CPP) -o $@ testrand.C

testrand$(EXT).out : testrand$(EXT)
	$? > $@


