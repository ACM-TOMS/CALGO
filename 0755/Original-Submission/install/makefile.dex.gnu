# 
#  --------------------------------------------------------------
#  makefile (DEX) of ADOL-C version 1.6 as of January 1,   1995
#  --------------------------------------------------------------
#  makefile for documented examples in subdirectory DEX 
#  written for GNU's  g++ compiler.
#
AD = ../
# AD may be any directory with ADOL-C library and header files 
#
CFLAG =  -O -I$(AD)
LFLAG = -L$(AD)  
CC = g++	
LIBS = -lad -lm 
all:    vectexam scalexam detexam odexam gaussexam

vectexam : vectexam.o $(AD)/libad.a
	$(CC) -o vectexam vectexam.o $(LFLAG) $(LIBS) 
vectexam.o : vectexam.c $(AD)/usrparms.h $(AD)/adouble.h $(AD)/adutils.h
	$(CC) -c $(CFLAG) vectexam.c

scalexam : scalexam.o $(AD)/libad.a
	$(CC) -o scalexam scalexam.o $(LFLAG) $(LIBS) 
scalexam.o : scalexam.c $(AD)/usrparms.h $(AD)/adouble.h $(AD)/adutils.h
	$(CC) -c $(CFLAG) scalexam.c

detexam : detexam.o $(AD)/libad.a
	$(CC) -g -o detexam detexam.o $(LFLAG) $(LIBS) 
detexam.o : detexam.c $(AD)/usrparms.h $(AD)/adouble.h $(AD)/adutils.h
	$(CC) -c $(CFLAG) detexam.c

odexam : odexam.o $(AD)/libad.a
	$(CC) -o odexam odexam.o $(LFLAG) $(LIBS) 
odexam.o : odexam.c $(AD)/usrparms.h $(AD)/adouble.h $(AD)/adutils.h
	$(CC) -c $(CFLAG) odexam.c

gaussexam : gaussexam.o $(AD)/libad.a
	$(CC) -o gaussexam gaussexam.o $(LFLAG) $(LIBS) 
gaussexam.o : gaussexam.c $(AD)/usrparms.h $(AD)/adouble.h $(AD)/adutils.h
	$(CC) -c $(CFLAG) gaussexam.c

clean:
	/bin/rm *.o
reallyclean:
	/bin/rm *exam
     
