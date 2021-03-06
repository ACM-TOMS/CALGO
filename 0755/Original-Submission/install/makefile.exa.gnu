#
#  --------------------------------------------------------------
#  makefile (EXA) of ADOL-C version 1.6 as of January 1,   1995
#  --------------------------------------------------------------
#  makefile for undocumented examples in subdirectory EXA
#  written for  GNU's  gcc and g++ compilers.
#
AD = ../
# AD may be any directory with ADOL-C library and header files 
#
CFLAG =  -O -I$(AD)
LFLAG = -L$(AD)   
CC = g++
MCC = gcc
LIBS = -lad -lm 
all:    vectexam scalexam detexam odexam gaussexam helm-diff-exam helm-auto-exam helm-vect-exam shuttlexam od2exam

helm-diff-exam  : helm-diff-exam.c
	$(MCC)  -o helm-diff-exam helm-diff-exam.c -lm

helm-auto-exam  : helm-auto-exam.o
	$(CC)  -o helm-auto-exam helm-auto-exam.o $(LFLAG) $(LIBS) 
helm-auto-exam.o: helm-auto-exam.c  $(AD)/adouble.h $(AD)/adutils.h
	$(CC) -c $(CFLAG) helm-auto-exam.c 

helm-vect-exam  : helm-vect-exam.o
	$(CC)  -o helm-vect-exam helm-vect-exam.o $(LFLAG) $(LIBS) 
helm-vect-exam.o: helm-vect-exam.c  $(AD)/adouble.h $(AD)/adutils.h
	$(CC) -c $(CFLAG) helm-vect-exam.c 

shuttlexam  : shuttlexam.o $(AD)/libad.a
	$(CC)  -o shuttlexam shuttlexam.o $(LFLAG) $(LIBS) 
shuttlexam.o: shuttlexam.c  $(AD)/usrparms.h $(AD)/adouble.h $(AD)/adutils.h
	$(CC) -c  $(CFLAG) shuttlexam.c 

vectexam : vectexam.o $(AD)/libad.a
	$(CC) -o vectexam vectexam.o $(LFLAG) $(LIBS) 
vectexam.o : vectexam.c $(AD)/usrparms.h $(AD)/adouble.h $(AD)/adutils.h
	$(CC) -c $(CFLAG) vectexam.c

scalexam : scalexam.o $(AD)/libad.a
	$(CC) -o scalexam scalexam.o $(LFLAG) $(LIBS) 
scalexam.o : scalexam.c $(AD)/usrparms.h $(AD)/adouble.h $(AD)/adutils.h
	$(CC) -c $(CFLAG) scalexam.c

detexam : detexam.o $(AD)/libad.a
	$(CC)  -o detexam detexam.o $(LFLAG) $(LIBS) 
detexam.o : detexam.c $(AD)/usrparms.h $(AD)/adouble.h $(AD)/adutils.h
	$(CC) -c $(CFLAG) detexam.c

odexam : odexam.o $(AD)/libad.a
	$(CC) -o odexam odexam.o $(LFLAG) $(LIBS) 
odexam.o : odexam.c $(AD)/usrparms.h $(AD)/adouble.h $(AD)/adutils.h
	$(CC) -c $(CFLAG) odexam.c

od2exam : od2exam.o eutroph.o $(AD)/libad.a
	$(CC) -o od2exam od2exam.o eutroph.o $(LFLAG) $(LIBS) 
od2exam.o : od2exam.c $(AD)/usrparms.h $(AD)/adouble.h $(AD)/adutils.h
	$(CC) -c $(CFLAG) od2exam.c
eutroph.o : eutroph.c $(AD)/usrparms.h $(AD)/adouble.h
	$(CC) -c $(CFLAG) eutroph.c

gaussexam : gaussexam.o $(AD)/libad.a
	$(CC) -o gaussexam gaussexam.o $(LFLAG) $(LIBS) 
gaussexam.o : gaussexam.c $(AD)/usrparms.h $(AD)/adouble.h $(AD)/adutils.h
	$(CC) -c $(CFLAG) gaussexam.c

clean:
	/bin/rm *.o
reallyclean:
	/bin/rm *exam
     
