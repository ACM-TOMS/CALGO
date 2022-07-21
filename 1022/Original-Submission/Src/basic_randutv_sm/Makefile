###############################################################################
#
# Makefile
#
###############################################################################

include $(INCLUDE_FILE)

#
# Define the C andlFortran90 compiler and loader, and their flags.
#
CC     = icc
FC     = gfortran
LINKER = icc

#
# Define the libraries.
#
LIBFLAME = /opt/libflame/

# General flags for compilation:
CFLAGS  += -O3 -qopenmp -I ${LIBFLAME}/include/ -DFLA_ENABLE_SUPERMATRIX
LDFLAGS += -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -pthread -lgfortran -static -lquadmath -O3 -qopenmp
LDLIBS  += ${LIBFLAME}/lib/libflame.a

#
# Define objects and dependencies.
#
UTV = \
      FLA_UTV_AB_UT_var49h.o \
      CPU_Inner_stubs_utv.o \
      MyFLA_Tools_QR_UT_var102.o \
      MyFLA_Normal_random.o \
      MyFLA_Parallel_normal_random.o \
      FLA_QR_Apply_QT.o \
      check_residuals.o \
      compute_svd.o \
      MyFLA_Utils.o \
      rank_sv.o \
      dlaputils.o \
      matrix_generate.o \
      dgemg.o \
      FLASH_Queue.o \
      flush_cache.o

BASIC_TEST = basic_test.o

all : basic_test.x

basic_test.x : $(BASIC_TEST) $(UTV) Makefile
	$(LINKER) $(BASIC_TEST) $(UTV) $(LDLIBS) $(LDFLAGS) \
                  -o basic_test.x

clean : 
	rm -f *.o *~ core *.x *__genmod.f90 *__genmod.mod
