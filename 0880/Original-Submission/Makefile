###########################################################################
#                                                                         #
# Usage:                                                                  #
#                                                                         #
#    make [target]                                                        #
#                                                                         #
# where target can be                                                     #
#                                                                         #
#    dstetester.x : double precision version of stetester                 #
#    sstetester.x : single precision version of stetester                 #
#    stetesterlib : creates the stetester library                         #
#                                                                         #
# By default all the three targets are created.                           #
#                                                                         #
###########################################################################

################
include make.inc
################

.SUFFIXES:
.SUFFIXES: .F .o

all: stetesterlib dstetester.x sstetester.x 

# Object files for stetester ##############################################

dstetester_lib = SRC/dstedefinitions.o \
		 SRC/dsteinterfaces1.o \
		 SRC/dsteinterfaces2.o \
		 SRC/dsteinterfaces3.o \
		 SRC/dstebndgap.o \
		 SRC/dstechkrslt.o \
		 SRC/dstedefeigv.o \
		 SRC/dstedefmtrx.o \
		 SRC/dstegetdata.o \
		 SRC/dstematgen.o \
		 SRC/dstemccallst.o \
		 SRC/dstemcdump.o \
		 SRC/dstemceigval.o \
		 SRC/dstemceigvalf.o \
		 SRC/dstemceigvi.o \
		 SRC/dstemceigvv.o \
		 SRC/dstemcglued.o \
	         SRC/dstemcmatrix.o \
		 SRC/dstemcmatrixf.o \
		 SRC/dsteprnrslt.o \
		 SRC/dsterndiliu.o \
		 SRC/dsterndvlvu.o \
		 SRC/dsteruntests.o \
		 SRC/getdreal.o

sstetester_lib = SRC/sstedefinitions.o \
		 SRC/ssteinterfaces1.o \
		 SRC/ssteinterfaces2.o \
		 SRC/ssteinterfaces3.o \
		 SRC/sstebndgap.o \
		 SRC/sstechkrslt.o \
		 SRC/sstedefeigv.o \
		 SRC/sstedefmtrx.o \
		 SRC/sstegetdata.o \
		 SRC/sstematgen.o \
		 SRC/sstemccallst.o \
		 SRC/sstemcdump.o \
		 SRC/sstemceigval.o \
		 SRC/sstemceigvalf.o \
		 SRC/sstemceigvi.o \
		 SRC/sstemceigvv.o \
		 SRC/sstemcglued.o \
		 SRC/sstemcmatrix.o \
		 SRC/sstemcmatrixf.o \
		 SRC/ssteprnrslt.o \
		 SRC/ssterndiliu.o \
		 SRC/ssterndvlvu.o \
		 SRC/ssteruntests.o \
		 SRC/getsreal.o

gstetester_lib = SRC/gstedefinitions.o \
		 SRC/gsteinterfaces1.o \
		 SRC/gsteinterfaces2.o \
		 SRC/getintgr.o \
		 SRC/getrecord.o \
		 SRC/getstrng.o \
		 SRC/handler.o \
		 SRC/listlenght.o \
		 SRC/parser.o \
		 SRC/parserlist.o

dstetester_obj = SRC/dstetester.o SRC/dsteflops.o SRC/dstetimer.o
sstetester_obj = SRC/sstetester.o SRC/ssteflops.o SRC/sstetimer.o

# Target _stetester.x requires LAPACK #####################################

dstetester.x: $(dstetester_obj)
	$(FC) $(FCLDFLGS) -o $@ $(dstetester_obj) $(STETESTERLIB) \
	$(PAPILIB) $(LAPACKLIB) $(TMGLIB) $(BLASLIB)

sstetester.x: $(sstetester_obj)
	$(FC) $(FCLDFLGS) -o $@ $(sstetester_obj) $(STETESTERLIB) \
	$(PAPILIB) $(LAPACKLIB) $(TMGLIB) $(BLASLIB)

# Target stetesterlib #####################################################

stetesterlib: $(gstetester_lib) $(dstetester_lib) $(sstetester_lib)
	$(ARCH) $(ARCHFLAGS) $(STETESTERLIB) \
		$(gstetester_lib) \
	        $(dstetester_lib) \
	        $(sstetester_lib)
	$(RANLIB) $(STETESTERLIB)

###########################################################################

clean:
	-@rm -f SRC/*.o *.log *.m *.out.* *.mod

cleanall:
	-@rm -f *.x LIB/*.a LIB/*.mod; make clean

.F.o:
	$(FC) $(FCFLAGS) $(FCMODLS) $(FCFRMT) -c -o $@ $<

SRC/dstetimer.o:
	$(FC) $(FCFLAGS) $(FCMODLS) $(FCFRMT) $(PAPIINC) \
              -c -o $@ SRC/dstetimer_$(TIMER).F
SRC/sstetimer.o:
	$(FC) $(FCFLAGS) $(FCMODLS) $(FCFRMT) $(PAPIINC) \
              -c -o $@ SRC/sstetimer_$(TIMER).F
SRC/dsteflops.o:
	$(FC) $(FCFLAGS) $(FCMODLS) $(FCFRMT) $(PAPIINC) \
              -c -o $@ SRC/dsteflops_$(FLOPS).F
SRC/ssteflops.o:
	$(FC) $(FCFLAGS) $(FCMODLS) $(FCFRMT) $(PAPIINC) \
              -c -o $@ SRC/ssteflops_$(FLOPS).F
