#######################################################################
#
# The user must set this options:
#

FORTRAN    = /opt/fortran/bin/f77
OPTS       = -u -g -C #-dalign
LOADER     = /opt/fortran/bin/f77
LOADOPTS   = -g -C +U77 #-dalign

LAPACKLIB  = -llapack-2 #/usr/local/lib/lapack/lapack.a
TMGLIB     = -ltmglib-2 #/usr/local/lib/lapack/tmglib.a
BLASLIB    = /opt/fortran/lib/pa1.1/libblas.a  #/usr/local/lib/lapack/blas.f77.a


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
#BLASLIB    = -lblas

#
#######################################################################


# No more changes are required beyond this line.


#######################################################################
#  Special modules
#  

RRQRLIB       = ../../rrqr.a
SHELL         = /bin/csh



ALINTST = \
   aladhd.o alaerh.o alaesm.o alahd.o alareq.o \
   alasum.o alasvm.o chkxer.o ilaenv.o xlaenv.o xerbla.o

SCLNTST= \
   slaord.o

DZLNTST= \
   dlaord.o 

SLINTST = schkaa.o schkqp.o schkrr.o \
   serrqp.o serrrr.o sqpt01.o sqrt11.o sqrt12.o srrt01.o srrt02.o

CLINTST = cchkaa.o cchkqp.o cchkrr.o \
   cerrqp.o cerrrr.o cqpt01.o cqrt11.o cqrt12.o crrt01.o crrt02.o

DLINTST = dchkaa.o dchkqp.o dchkrr.o \
   derrqp.o derrrr.o dqpt01.o dqrt11.o dqrt12.o drrt01.o drrt02.o

ZLINTST = zchkaa.o zchkqp.o zchkrr.o \
   zerrqp.o zerrrr.o zqpt01.o zqrt11.o zqrt12.o zrrt01.o zrrt02.o

all:  single double complex complex16

single: ../xlintsts
double: ../xlintstd
complex: ../xlintstc
complex16: ../xlintstz

../xlintsts : $(ALINTST) $(SLINTST) $(SCLNTST)
	$(LOADER) $(LOADOPTS)  $(ALINTST) $(SCLNTST) $(SLINTST) \
                  $(RRQRLIB) $(TMGLIB) $(LAPACKLIB) $(BLASLIB)  -o $@

../xlintstc : $(ALINTST) $(CLINTST) $(SCLNTST)
	$(LOADER) $(LOADOPTS)  $(ALINTST) $(SCLNTST) $(CLINTST) \
                  $(RRQRLIB) $(TMGLIB) $(LAPACKLIB) $(BLASLIB)  -o $@
 
../xlintstd : $(ALINTST) $(DLINTST) $(DZLNTST)
	$(LOADER) $(LOADOPTS)  $(ALINTST) $(DZLNTST) $(DLINTST) \
                  $(RRQRLIB) $(TMGLIB) $(LAPACKLIB) $(BLASLIB)  -o $@

../xlintstz : $(ALINTST) $(ZLINTST) $(DZLNTST)
	$(LOADER) $(LOADOPTS)  $(ALINTST) $(DZLNTST) $(ZLINTST) \
                  $(RRQRLIB) $(TMGLIB) $(LAPACKLIB) $(BLASLIB)  -o $@

#
# Rules for the modules.
#

.f.o: 
	$(FORTRAN) -c $(OPTS) $*.f

#
# Clean rule.
#

clean:
	- rm -f *.o ../xlintst?
