#######################################################################
#
# The user must set this options:
#

FORTRAN   = f77
OPTS      = -u -g -C

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

#
# On SUN computers, it is recommended:
#
#FORTRAN   = f77
#OPTS      = -u -O
#
# On IBM RS6K, it is recommended:
#
#FORTRAN   = xlf
#OPTS      = -u -O3
#

#
#######################################################################


# No more changes are required beyond this line.


#######################################################################
#
# Other definitions:
#

REAL_SOURCES = ../v15.acm
COMPLEX_SOURCES = ../cv15.acm

RRQR_LIB = ../rrqr.a

#
#######################################################################

#######################################################################
#  
#  Modules for Rank-Revealing QR:
#  

S_RRQR_MODULES = \
   sgeqpb.o sgeqpw.o sgeqpc.o \
   sgeqpx.o strqpx.o strqxc.o \
   sgeqpy.o strqpy.o strqyc.o \
   strrnk.o slauc1.o slasmx.o \
   smylap.o

D_RRQR_MODULES = \
   dgeqpb.o dgeqpw.o dgeqpc.o \
   dgeqpx.o dtrqpx.o dtrqxc.o \
   dgeqpy.o dtrqpy.o dtrqyc.o \
   dtrrnk.o dlauc1.o dlasmx.o \
   dmylap.o

C_RRQR_MODULES = \
   cgeqpb.o cgeqpw.o cgeqpc.o \
   cgeqpx.o ctrqpx.o ctrqxc.o \
   cgeqpy.o ctrqpy.o ctrqyc.o \
   ctrrnk.o clauc1.o clasmx.o \
   cmylap.o

Z_RRQR_MODULES = \
   zgeqpb.o zgeqpw.o zgeqpc.o \
   zgeqpx.o ztrqpx.o ztrqxc.o \
   zgeqpy.o ztrqpy.o ztrqyc.o \
   ztrrnk.o zlauc1.o zlasmx.o \
   zmylap.o

#  
#  
#######################################################################


all:  single double complex complex16

single: $(S_RRQR_MODULES)
	$(ARCH) $(ARCHFLAGS) $(RRQR_LIB) $(S_RRQR_MODULES)
	$(RANLIB) $(RRQR_LIB)

double: $(D_RRQR_MODULES)
	$(ARCH) $(ARCHFLAGS) $(RRQR_LIB) $(D_RRQR_MODULES)
	$(RANLIB) $(RRQR_LIB)

complex: $(C_RRQR_MODULES)
	$(ARCH) $(ARCHFLAGS) $(RRQR_LIB) $(C_RRQR_MODULES)
	$(RANLIB) $(RRQR_LIB)

complex16: $(Z_RRQR_MODULES)
	$(ARCH) $(ARCHFLAGS) $(RRQR_LIB) $(Z_RRQR_MODULES)
	$(RANLIB) $(RRQR_LIB)

clean:
	- rm -f *.o

#
# Rules for the object files of RRQR.
#

.f.o:

#
# Rules for the source modules of RRQR.
#

sgeqpb.o: $(REAL_SOURCES)/sgeqpb.f
	$(FORTRAN) -c $(OPTS) $?

sgeqpw.o: $(REAL_SOURCES)/sgeqpw.f
	$(FORTRAN) -c $(OPTS) $?

sgeqpc.o: $(REAL_SOURCES)/sgeqpc.f
	$(FORTRAN) -c $(OPTS) $?

sgeqpx.o: $(REAL_SOURCES)/sgeqpx.f
	$(FORTRAN) -c $(OPTS) $?

strqpx.o: $(REAL_SOURCES)/strqpx.f
	$(FORTRAN) -c $(OPTS) $?

strqxc.o: $(REAL_SOURCES)/strqxc.f
	$(FORTRAN) -c $(OPTS) $?

sgeqpy.o: $(REAL_SOURCES)/sgeqpy.f
	$(FORTRAN) -c $(OPTS) $?

strqpy.o: $(REAL_SOURCES)/strqpy.f
	$(FORTRAN) -c $(OPTS) $?

strqyc.o: $(REAL_SOURCES)/strqyc.f
	$(FORTRAN) -c $(OPTS) $?

strrnk.o: $(REAL_SOURCES)/strrnk.f
	$(FORTRAN) -c $(OPTS) $?

slauc1.o: $(REAL_SOURCES)/slauc1.f
	$(FORTRAN) -c $(OPTS) $?

slasmx.o: $(REAL_SOURCES)/slasmx.f
	$(FORTRAN) -c $(OPTS) $?

smylap.o: $(REAL_SOURCES)/smylap.f
	$(FORTRAN) -c $(OPTS) $?


dgeqpb.o: $(REAL_SOURCES)/dgeqpb.f
	$(FORTRAN) -c $(OPTS) $?

dgeqpw.o: $(REAL_SOURCES)/dgeqpw.f
	$(FORTRAN) -c $(OPTS) $?

dgeqpc.o: $(REAL_SOURCES)/dgeqpc.f
	$(FORTRAN) -c $(OPTS) $?

dgeqpx.o: $(REAL_SOURCES)/dgeqpx.f
	$(FORTRAN) -c $(OPTS) $?

dtrqpx.o: $(REAL_SOURCES)/dtrqpx.f
	$(FORTRAN) -c $(OPTS) $?

dtrqxc.o: $(REAL_SOURCES)/dtrqxc.f
	$(FORTRAN) -c $(OPTS) $?

dgeqpy.o: $(REAL_SOURCES)/dgeqpy.f
	$(FORTRAN) -c $(OPTS) $?

dtrqpy.o: $(REAL_SOURCES)/dtrqpy.f
	$(FORTRAN) -c $(OPTS) $?

dtrqyc.o: $(REAL_SOURCES)/dtrqyc.f
	$(FORTRAN) -c $(OPTS) $?

dtrrnk.o: $(REAL_SOURCES)/dtrrnk.f
	$(FORTRAN) -c $(OPTS) $?

dlauc1.o: $(REAL_SOURCES)/dlauc1.f
	$(FORTRAN) -c $(OPTS) $?

dlasmx.o: $(REAL_SOURCES)/dlasmx.f
	$(FORTRAN) -c $(OPTS) $?

dmylap.o: $(REAL_SOURCES)/dmylap.f
	$(FORTRAN) -c $(OPTS) $?


cgeqpb.o: $(COMPLEX_SOURCES)/cgeqpb.f
	$(FORTRAN) -c $(OPTS) $?

cgeqpw.o: $(COMPLEX_SOURCES)/cgeqpw.f
	$(FORTRAN) -c $(OPTS) $?

cgeqpc.o: $(COMPLEX_SOURCES)/cgeqpc.f
	$(FORTRAN) -c $(OPTS) $?

cgeqpx.o: $(COMPLEX_SOURCES)/cgeqpx.f
	$(FORTRAN) -c $(OPTS) $?

ctrqpx.o: $(COMPLEX_SOURCES)/ctrqpx.f
	$(FORTRAN) -c $(OPTS) $?

ctrqxc.o: $(COMPLEX_SOURCES)/ctrqxc.f
	$(FORTRAN) -c $(OPTS) $?

cgeqpy.o: $(COMPLEX_SOURCES)/cgeqpy.f
	$(FORTRAN) -c $(OPTS) $?

ctrqpy.o: $(COMPLEX_SOURCES)/ctrqpy.f
	$(FORTRAN) -c $(OPTS) $?

ctrqyc.o: $(COMPLEX_SOURCES)/ctrqyc.f
	$(FORTRAN) -c $(OPTS) $?

ctrrnk.o: $(COMPLEX_SOURCES)/ctrrnk.f
	$(FORTRAN) -c $(OPTS) $?

clauc1.o: $(COMPLEX_SOURCES)/clauc1.f
	$(FORTRAN) -c $(OPTS) $?

clasmx.o: $(COMPLEX_SOURCES)/clasmx.f
	$(FORTRAN) -c $(OPTS) $?

cmylap.o: $(COMPLEX_SOURCES)/cmylap.f
	$(FORTRAN) -c $(OPTS) $?


zgeqpb.o: $(COMPLEX_SOURCES)/zgeqpb.f
	$(FORTRAN) -c $(OPTS) $?

zgeqpw.o: $(COMPLEX_SOURCES)/zgeqpw.f
	$(FORTRAN) -c $(OPTS) $?

zgeqpc.o: $(COMPLEX_SOURCES)/zgeqpc.f
	$(FORTRAN) -c $(OPTS) $?

zgeqpx.o: $(COMPLEX_SOURCES)/zgeqpx.f
	$(FORTRAN) -c $(OPTS) $?

ztrqpx.o: $(COMPLEX_SOURCES)/ztrqpx.f
	$(FORTRAN) -c $(OPTS) $?

ztrqxc.o: $(COMPLEX_SOURCES)/ztrqxc.f
	$(FORTRAN) -c $(OPTS) $?

zgeqpy.o: $(COMPLEX_SOURCES)/zgeqpy.f
	$(FORTRAN) -c $(OPTS) $?

ztrqpy.o: $(COMPLEX_SOURCES)/ztrqpy.f
	$(FORTRAN) -c $(OPTS) $?

ztrqyc.o: $(COMPLEX_SOURCES)/ztrqyc.f
	$(FORTRAN) -c $(OPTS) $?

ztrrnk.o: $(COMPLEX_SOURCES)/ztrrnk.f
	$(FORTRAN) -c $(OPTS) $?

zlauc1.o: $(COMPLEX_SOURCES)/zlauc1.f
	$(FORTRAN) -c $(OPTS) $?

zlasmx.o: $(COMPLEX_SOURCES)/zlasmx.f
	$(FORTRAN) -c $(OPTS) $?

zmylap.o: $(COMPLEX_SOURCES)/zmylap.f
	$(FORTRAN) -c $(OPTS) $?
