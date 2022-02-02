
#************************************************************************
# Author:	Leonidas Linardakis
#
#
# File:		Makefile
# Date:		03 Sep 2005
# Description: creates the madd 	
#
#************************************************************************

include ../make.settings
CC=${CompCplusplus}

OBJS = madd.o maddi.o L_Timer.o


default: all
all: maddi libweight 

maddi: $(OBJS) Makefile 
	${CC} ${Optflags} ${OBJS}  ${linktriangle}  ${linkmetis} -ldl -o $@   

maddi.o: maddi.cc madd.h 
	${CC} ${Optflags} -c $< 

madd.o: madd.cc madd.h 
	${CC} ${Optflags} -c $<

L_Timer.o: L_Timer.cc L_Timer.h
	${CC} ${Optflags} -c  $<

libweight: weightLibrary.o
	${dlCC} -shared  weightLibrary.o  -o libweight.so -lc

weightLibrary.o: weightLibrary.cc
	${dlCC} -fPIC -c  $<







