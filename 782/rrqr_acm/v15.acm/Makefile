#
#      Makefile for generating single and double precision routines
#      for rrqr factorization.
#
#      ***********************************************
#      * Authors:  C.H.Bischof  and  G.Quintana-Orti *
#      *                                             *
#      ***********************************************
#
# To be set in command line:
# =========================

.SUFFIX:

P         = s                    # Precision (s or d).
FORTRAN   = f77	                 # Fortran compiler.
OPTS      = -u -dalign           # Fortran compiler optimizing options.
LOADER    = f77                  # Fortran linker.
LOADOPTS  = -O -dalign           # Fortran linker optimizing options.
LAPACKLIB = /usr/local/lib/lapack/lapack-2.a    #  Name of libraries 
						# supplying LAPACK 
						# solvers.
TMGLIB    = /usr/local/lib/lapack/tmglib-2.a   	#  Name of libraries 
						# supplying LAPACK 
						# test matrix 
						# generation suite.
BLASLIB   = /usr/local/lib/lapack/blas.a        #  Name of libraries 
						# supplying BLAS.
 
SHELL = /bin/csh

OBJECTS = $(P)qr.o \
	$(P)geqpx.o $(P)trqpx.o $(P)trqxc.o \
	$(P)geqpy.o $(P)trqpy.o $(P)trqyc.o \
	$(P)geqpb.o $(P)geqpw.o $(P)geqpc.o \
	$(P)mylap.o \
	$(P)trrnk.o $(P)lasmx.o $(P)lauc1.o \
	$(P)gntst.o $(P)qrmtx.o \
	$(P)utils.o \
	ilaenv.o \
	esm.o

#      -- Main Target --
#      ******************

$(P)qr: $(OBJECTS)
	$(LOADER) $(LOADOPTS) -o $@ $(OBJECTS) \
                                    $(TMGLIB) $(LAPACKLIB) $(BLASLIB)

#      -- Rules for generating objects --
#      **********************************

qr.o:
	echo $(P)qr

.f.o:
	$(FORTRAN) -c $(OPTS) $*.f

#      -- Rules for cleaning the source and object files --
#      ****************************************************

clean: cleanup cleanout

cleanup:
	- rm -f sqr_* dqr_* *.o

cleanout:
	- rm -f out out.* time.* rank.* ftrace.* *.trace core
