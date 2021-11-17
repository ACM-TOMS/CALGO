#-------------------------------------------------------------------------------
# Makefile for the LDL mexFunction
#-------------------------------------------------------------------------------

# edit the following three lines to select your compiler and archiver
CC = cc

# with picky warnings and statement coverage:
# CC = gcc -Wall -W -Wshadow -Wmissing-prototypes -Wstrict-prototypes \
	-pedantic -ansi -g -ftest-coverage -fprofile-arcs

AR = ar
RANLIB = ranlib

all: intro libldl.a ldlsimple ldlmain ldl_mexfunction ldlamd

intro:
	@echo " "
	@echo "    Now compiling the LDL library (libldl.a), and three stand-"
	@echo "    alone C programs (ldlmain, ldlsimple, and ldlamd).  ldlamd"
	@echo "    is ldlmain compiled to use with AMD, which will be skipped"
	@echo "    if you do not have AMD installed in the ../AMD directory."
	@echo "    Will also compile the LDL mexFunction for use in MATLAB,"
	@echo "    and two MATLAB test programs, which will be skipped if you"
	@echo "    do not have MATLAB."
	@echo " "

#-------------------------------------------------------------------------------
# the ldl library:
#-------------------------------------------------------------------------------

libldl.a: ldl.c ldl.h
	$(CC) -c -O ldl.c -o ldl.o
	$(AR) cr libldl.a ldl.o
	- $(RANLIB) libldl.a

#-------------------------------------------------------------------------------
# stand-alone C programs:
#-------------------------------------------------------------------------------

ldlmain:  ldlmain.c libldl.a
	$(CC) -O ldlmain.c libldl.a -o ldlmain -lm
	- ldlmain > my_ldlmain.out
	- diff ldlmain.out my_ldlmain.out

ldlsimple:  ldlsimple.c libldl.a
	$(CC) -O ldlsimple.c libldl.a -o ldlsimple -lm
	- ldlsimple > my_ldlsimple.out
	- diff ldlsimple.out my_ldlsimple.out

ldlamd:  ldlmain.c libldl.a
	- make -C ../AMD
	- $(CC) -O -I../AMD/Include -L../AMD/Lib -DUSE_AMD ldlmain.c -lamd libldl.a -o ldlamd -lm
	- ldlamd > my_ldlamd.out
	- diff ldlamd.out my_ldlamd.out

#-------------------------------------------------------------------------------
# LDL mexFunctions for use in MATLAB:
#-------------------------------------------------------------------------------

ldl_mexfunction: ldlmex.c ldl.c ldl.h ldlsymbolmex.c
	- mex -inline ldl.c ldlmex.c
	- mex -inline -output ldlsymbol ldlsymbolmex.c ldl.c

#-------------------------------------------------------------------------------
# clean-up:
#-------------------------------------------------------------------------------

clean:
	- 'rm' -f *.o

distclean: purge

purge: clean
	- 'rm' ldl.mex* ldl.dll libldl.a ldlmain ldlamd ldlsimple
	- 'rm' my_ldlmain.out my_ldlamd.out my_ldlsimple.out
	- 'rm' ldlmain.mex* ldlamd.mex* ldlmain.dll ldlamd.dll
	- 'rm' ldlsymbol.mex* ldlsymbol.dll
	- 'rm' *.dvi *.aux *.log *.bak *.bbl *.blg


#-------------------------------------------------------------------------------
# user guide:
#-------------------------------------------------------------------------------

doc: ldl_userguide.tex ldl.bib
	latex ldl_userguide
	- bibtex ldl_userguide
	latex ldl_userguide
	latex ldl_userguide
	dvips ldl_userguide -o ldl_userguide.ps
	pdflatex ldl_userguide
	pdflatex ldl_userguide

