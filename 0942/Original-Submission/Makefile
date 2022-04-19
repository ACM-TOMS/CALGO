# Makefile for StencilProbe
#
# Modified by Raúl de la Cruz to add Semi-stencil and hwc versions support
# Computer Application department (CASE)
# Barcelona Supercomputing Center (BSC) - Spain
#
# Valid parameters for this Makefile:
# OPTS: Can contain those flags that must be passed to the compiler refering
#       to source code modifications/optimizations. Valid modifiers are:
#       -DFISSION_2LOOPS, -DFISSION_3LOOPS: apply fission to internal loop
#       -DSEMI_ALL: apply Semi-stencil computation algorithm into three axis
#       -DNUM_TRIALS: change the number of trials to get performance results
#       TRACE: tracing support for HWC binaries
#       	-DSEQTRACE: use EXTRAE package (useful in BG/P)
#       	-DPAPITRACE: use PAPI library directly, set PAPI_COUNTERS environment
#       	             variable with PAPI presets in a comma-separated fashion
#       	NOTHING: use PAPIEX tool
#       -DPLOT: pretty printing when getting performance results (times & hwc)
#       -DPEAK: remove initializations and clearing cache in order to measure
#               actual peak performance
#
# SUFFIX_DEF: suffix name for default binaries (to get performance results)
# SUFFIX_HWC: suffix name for hardware counters binaries (to get metrics)
#
# OPENMP: enable openmp implementation of the stencilprobe algorithms
#

###################################################################
# Select architecture (AMD,XEON,POWER6,POWER7,BGP,SANDY,MIC) ######
ARCH=MIC
###################################################################

# Compilers and flags used on every architecture
#

#
# Sandy Bridge: -DAVX macro is enabled by default (vectorial code),
# remove it to use general purpose code (ICC preferred)
CC_SANDY      = icc
CFLAGS_SANDY  = -O3 -xHOST -ipo -no-prec-div -vec-report=6 -fno-alias -g $(OPTS) -DAVX
OMPFLAG_SANDY = -openmp
#CC_SANDY      = gcc
#CFLAGS_SANDY  = -O3 -g $(OPTS)
#OMPFLAG_SANDY = -fopenmp

#
# Louhi: AMD64 Quad-Core AMD Opteron(tm) Processor 23 (C2) 2.7GHz (prace-shangai queue) (GCC preferred)
CC_AMD      = gcc
CFLAGS_AMD  = -O3 -march=amdfam10 -mtune=amdfam10 -g $(OPTS)
OMPFLAG_AMD = -fopenmp
# CC_AMD      = pgcc
# CFLAGS_AMD  = -O3 -tp shanghai-64 -g $(OPTS)
# OMPFLAG_AMD = -mp

#
# Inti: x86_64 Intel(R) Xeon(R) CPU X5570 @ 2.93GHz (GCC preferred)
# Juropa: x86_64 Intel(R) Xeon(R) CPU X5570 @ 2.93GHz (GCC preferred)
CC_XEON      = gcc
CFLAGS_XEON  = -O3 -mtune=core2 -march=core2 -msse4.1 -g $(OPTS)
OMPFLAG_XEON = -fopenmp
# CC_XEON      = icc
# CFLAGS_XEON  = -O3 -mtune=core2 -march=core2 -msse4.1 -g $(OPTS)
# OMPFLAG_XEON = -openmp

#
# Huygens: POWER6 (architected), altivec supported 4.7GHz (XLC preferred)
CC_PWR6      = xlc
CFLAGS_PWR6  = -O3 -qarch=pwr6 -qtune=pwr6 -g -I/sara/sw/papi/3.6.2/include $(OPTS)
OMPFLAG_PWR6 = -qsmp=omp -qthreaded

#
# ps701n1: POWER7 (architected), altivec supported 3.0GHz (XLC preferred)
CC_PWR7      = xlc
CFLAGS_PWR7  = -O3 -q64 -qarch=pwr7 -qtune=pwr7 -g -I/home/delacruz/local/include $(OPTS)
OMPFLAG_PWR7 = -qsmp=omp -qthreaded

#
# Jugene: BlueGene/P PowerPC 440 850Mhz (XL Compilers and GCC for oblivious bins)
CC_BGP      = mpixlc # bgcc
CFLAGS_BGP  = -O3 -g $(OPTS) -I~/local/include
OMPFLAG_BGP = -qsmp=omp -qthreaded
# CC_BGP      = mpicc # /bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc-bgp-linux-gcc
# CFLAGS_BGP  = -O3 -g $(OPTS) -I~/local/include
# OMPFLAG_BGP = -qsmp=omp -qthreaded

#
# bscin02 laptop i686: Intel(R) Core(TM)2 Duo CPU T9400 @ 2.53GHz (GCC preferred)
CC_BSCIN02      = gcc
CFLAGS_BSCIN02  = -O3 -mtune=core2 -march=core2 -g $(OPTS)
OMPFLAG_BSCIN02 = -fopenmp

#
# Knights Ferry/Corner: Intel(R) MIC 1GHz. -DMIC macro is enabled by default
# (vectorial code), remove it to use general purpose code (ICC preferred)
CC_MIC      = icc
CFLAGS_MIC  = -O3 -mmic -fno-alias -fargument-noalias -vec-report=3 -g $(OPTS) -DMIC # -opt-prefetch=0
OMPFLAG_MIC = -openmp


# Global flags
CC = $(CC_$(ARCH))
CFLAGS = $(CFLAGS_$(ARCH))
OMPFLAG = $(OMPFLAG_$(ARCH))


# Uncomment this line in order to use PAPI for HWC monitoring
#PAPI = -lpapi
ifndef SUFFIX_HWC
  PAPI =
endif
LDFLAGS = $(PAPI) #-lm

# the line below defines timers.  if not defined, will attempt to automatically
# detect available timers.  See cycle.h.
# should be set to -DHAVE_PAPI or -DHAVE_GETTIMEOFDAY or unset.
#TIMER = -DHAVE_PAPI
TIMER = -DHAVE_GETTIMEOFDAY


ifdef OPENMP
  CFLAGS+=$(OMPFLAG)
endif
CFLAGS += -I$(PWD)/includes


# Define ending name for benchmarking tests
ifndef SUFFIX_DEF
  SUFFIX_DEF=
endif

ifndef SUFFIX_HWC
  SUFFIX_HWC=.hwc
endif

VPATH = sources:includes

BIN = $(PWD)/bin

DEF = probe$(SUFFIX_DEF)      blocked_probe$(SUFFIX_DEF)      timeskew_probe$(SUFFIX_DEF)      oblivious_probe$(SUFFIX_DEF)      \
      semi_probe$(SUFFIX_DEF) blocked_semi_probe$(SUFFIX_DEF) timeskew_semi_probe$(SUFFIX_DEF) oblivious_semi_probe$(SUFFIX_DEF)

HWC = probe$(SUFFIX_HWC)      blocked_probe$(SUFFIX_HWC)      timeskew_probe$(SUFFIX_HWC)      oblivious_probe$(SUFFIX_HWC)      \
      semi_probe$(SUFFIX_HWC) blocked_semi_probe$(SUFFIX_HWC) timeskew_semi_probe$(SUFFIX_HWC) oblivious_semi_probe$(SUFFIX_HWC)

TEST = test$(SUFFIX_DEF)

BINS = $(DEF) $(HWC) $(TEST)

all: $(DEF)

def: $(DEF)

hwc: $(HWC)


#
# Default versions (with naive stencil computation)
#
probe$(SUFFIX_DEF): main.c util.c run.h probe_heat.c cycle.h stencil.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

blocked_probe$(SUFFIX_DEF): main.c util.c probe_heat_blocked.c cycle.h stencil.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

timeskew_probe$(SUFFIX_DEF): main.c util.c run.h probe_heat_timeskew.c cycle.h stencil.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES -DTIMESKEW_BLK $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

oblivious_probe$(SUFFIX_DEF): main.c util.c run.h probe_heat_oblivious.c cycle.h stencil.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

circqueue_probe$(SUFFIX_DEF): main.c util.c run.h probe_heat_circqueue.c cycle.h stencil.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES -DCIRCULARQUEUEPROBE $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

#
# Optimized versions (with semi-stencil computation)
#
semi_probe$(SUFFIX_DEF): main.c util.c probe_heat_semi.c cycle.h semi.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

blocked_semi_probe$(SUFFIX_DEF): main.c util.c probe_heat_blocked_semi.c cycle.h semi.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

timeskew_semi_probe$(SUFFIX_DEF): main.c util.c probe_heat_timeskew_semi.c cycle.h semi.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES -DTIMESKEW_BLK $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

oblivious_semi_probe$(SUFFIX_DEF): main.c util.c probe_heat_oblivious_semi.c cycle.h semi.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

#
# HWC default versions (used to get hardware counters)
#
probe$(SUFFIX_HWC): main.hwc.c util.c run.h probe_heat.c cycle.h stencil.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES -DHWC $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

blocked_probe$(SUFFIX_HWC): main.hwc.c util.c probe_heat_blocked.c cycle.h stencil.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES -DHWC $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

timeskew_probe$(SUFFIX_HWC): main.hwc.c util.c run.h probe_heat_timeskew.c cycle.h stencil.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES -DTIMESKEW_BLK -DHWC $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

oblivious_probe$(SUFFIX_HWC): main.hwc.c util.c run.h probe_heat_oblivious.c cycle.h stencil.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES -DHWC $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

#
# HWC optimized versions (used to get hardware counters)
#
semi_probe$(SUFFIX_HWC): main.hwc.c util.c probe_heat_semi.c cycle.h semi.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES -DHWC $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

blocked_semi_probe$(SUFFIX_HWC): main.hwc.c util.c probe_heat_blocked_semi.c cycle.h semi.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES -DHWC $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

timeskew_semi_probe$(SUFFIX_HWC): main.hwc.c util.c probe_heat_timeskew_semi.c cycle.h semi.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES -DTIMESKEW_BLK -DHWC $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

oblivious_semi_probe$(SUFFIX_HWC): main.hwc.c util.c probe_heat_oblivious_semi.c cycle.h semi.h
	$(CC) $(CFLAGS) $(TIMER) -DRANDOMVALUES -DHWC $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

#
# Testing program (checks if the computation is correct for each algorithm)
#
test$(SUFFIX_DEF): main.test.c util.c run.h probe_heat.c cycle.h probe_heat_blocked.c \
      probe_heat_oblivious.c probe_heat_timeskew.c probe_heat_circqueue.c     \
      probe_heat_semi.c probe_heat_blocked_semi.c probe_heat_oblivious_semi.c \
      probe_heat_timeskew_semi.c semi.h stencil.h
	$(CC) $(CFLAGS) $(TIMER) -DSTENCILTEST $(filter %.c,$^) $(LDFLAGS) -o $(BIN)/$@

clean:
	rm -f *.o $(BIN)/*

