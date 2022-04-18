#######################################################################
#
# The user must set this options:
#

FORTRAN    = f77
OPTS       = -u -O
LOADER     = f77 +U77
LOADOPTS   = -O

LAPACKLIB  = /usr/local/lib/lapack/lapack-2.a
TMGLIB     = /usr/local/lib/lapack/tmglib-2.a
BLASLIB    = /usr/local/lib/lapack/blas.f77.a


# On the SUN Solaris, it is recommended:
#FORTRAN    = f77
#OPTS       = -u -O -dalign
#LOADER     = f77
#LOADOPTS   = -O -dalign

#LAPACKLIB  =
#TMGLIB     = /usr/local/lapack/tmglib.a
#BLASLIB    = /home1/SUNWspro/SC3.0.1/lib/libsunperf.a


# On the IBM RS6K, it is recommended:
#FORTRAN    = xlf
#OPTS       = -u -O3
#LOADER     = xlf
#LOADOPTS   = -O3 -bnso -bI:/lib/syscalls.exp
#
#LAPACKLIB  = /usr/local/lapack/lib/lapack.a
#TMGLIB     = /usr/local/lapack/lib/tmglib.a
#BLASLIB    = -lessl -lblas

#
#######################################################################


# No more changes are required beyond this line.


#######################################################################
#  Special modules
#  

RRQRLIB       = ../../rrqr.a
SHELL         = /bin/csh

ALINTIM =  atimin.o atimck.o icopy.o ilaenv.o xlaenv.o

SCINTIM =  sprtbl.o sprtb4.o sprtb5.o \
         smflop.o sopbl2.o sopbl3.o sopla.o
 
DZINTIM =  dprtbl.o dprtb4.o dprtb5.o \
         dmflop.o dopbl2.o dopbl3.o dopla.o
 
SLINTIM = stimaa.o stimmg.o \
         stimmv.o stimmm.o stimqr.o stimqp.o stimrr.o
 
CLINTIM = ctimaa.o ctimmg.o \
         ctimmv.o ctimmm.o ctimqr.o ctimqp.o ctimrr.o
 
DLINTIM = dtimaa.o dtimmg.o \
         dtimmv.o dtimmm.o dtimqr.o dtimqp.o dtimrr.o
 
ZLINTIM = ztimaa.o ztimmg.o \
         ztimmv.o ztimmm.o ztimqr.o ztimqp.o ztimrr.o
 
all:  single double complex complex16

single: ../xlintims
double: ../xlintimd
complex: ../xlintimc
complex16: ../xlintimz

../xlintims : $(ALINTIM) $(SCINTIM) $(SLINTIM)
	$(LOADER) $(LOADOPTS)  $(ALINTIM) $(SCINTIM) $(SLINTIM) \
	          $(RRQRLIB) $(TMGLIB) $(LAPACKLIB) $(BLASLIB)  -o $@

../xlintimc : $(ALINTIM) $(SCINTIM) $(CLINTIM)
	$(LOADER) $(LOADOPTS)  $(ALINTIM) $(SCINTIM) $(CLINTIM) \
	          $(RRQRLIB) $(TMGLIB) $(LAPACKLIB) $(BLASLIB) -o $@

../xlintimd : $(ALINTIM) $(DZINTIM) $(DLINTIM)
	$(LOADER) $(LOADOPTS)  $(ALINTIM) $(DZINTIM) $(DLINTIM) \
	          $(RRQRLIB) $(TMGLIB) $(LAPACKLIB) $(BLASLIB)  -o $@

../xlintimz : $(ALINTIM) $(DZINTIM) $(ZLINTIM)
	$(LOADER) $(LOADOPTS)  $(ALINTIM) $(DZINTIM) $(ZLINTIM) \
	          $(RRQRLIB) $(TMGLIB) $(LAPACKLIB) $(BLASLIB)  -o $@

$(ALINTIM): $(FRC)
$(SCINTIM): $(FRC)
$(DZINTIM): $(FRC)
$(SLINTIM): $(FRC)
$(CLINTIM): $(FRC)
$(DLINTIM): $(FRC)
$(ZLINTIM): $(FRC)
 
FRC:
	@FRC=$(FRC)

 
#
# Rules for the modules.
#
.f.o : 
	$(FORTRAN) $(OPTS) -c $*.f

#
# Clean rule.
#

clean:
	- rm -f *.o ../xlintim?
