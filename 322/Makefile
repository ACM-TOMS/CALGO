#Makefile for fisher.c and driver program.
#Copyright (c) 1995 by Matthew Belmonte <mkb4@Cornell.edu>.
#Permission for use and distribution is hereby granted, subject to the
#restrictions that this copyright notice be included in its entirety, and that
#any and all changes are clearly noted.
#
CC = /bin/cc
CFLAGS = -O2

ftest: ftest.o fisher.o
	$(CC) $(CFLAGS) -o ftest ftest.o fisher.o -lm

ftest.o: ftest.c fisher.h
	$(CC) $(CFLAGS) -c ftest.c

fisher.o: fisher.c
	$(CC) $(CFLAGS) -c fisher.c

clean:
	rm -f ftest.o fisher.o ftest
