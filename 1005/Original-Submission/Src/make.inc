# -*- mode:makefile-*-
# See comments in Makefile in the blas-rmd root directory.
SHELL=/bin/bash
ifndef MAKE_INC
MAKE_INC:=true # prevent multiple inclusion

OS := $(shell uname -s)
ISMAC := $(filter $(OS),Darwin)
ifeq ($(origin FC),default)
  FC := gfortran
endif

ifeq ($(BLAS),) # $BLAS not set
  ifeq ($(blas),)
    BLAS := openblas
  else
    BLAS := $(blas)
  endif
endif
ifeq ($(OPENBLAS),)# $OPENBLAS not set
  OPENBLAS:=/usr/local/opt/openblas
endif

# See the file readme for explanation of the following:
gfortranmkl       := -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
gfortranstaticmkl := -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
ifortmkl          := -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
macmkl            := -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

# Blas-settings:
BLASDEP :=
ifeq ($(BLAS),openblas)
  BLASLIB := -L$(OPENBLAS)/lib -lopenblas
  NUMTHRD :=  num_threads_openblas
else ifeq ($(BLAS),mkl)
  ifeq ($(FC),gfortran)
    BLASLIB := $(if $(ISMAC),$(macmkl),$(gfortranmkl))
    FOPTS   += -m64
  else ifeq ($(FC),nagfor)
    BLASLIB := $(if $(ISMAC),$(macmkl),$(gfortranmkl))
  else ifeq ($(FC),ifort)
    BLASLIB := $(ifortmkl)
  else
    $(error bad Fortran)
  endif
  FOPTS += -I${MKLROOT}/include
  NUMTHRD :=  num_threads_mkl
else ifeq ($(BLAS),accelerate)
  BLASLIB := -framework accelerate
  NUMTHRD :=  num_threads_noop
else ifeq ($(BLAS),refblas)
  NUMTHRD := num_threads_noop
  ifeq ($(TOPLEVEL),true)
    BLASLIB := -Llib -lrefblas
    BLASDEP := lib/librefblas.a
  else
    BLASLIB := -L../lib -lrefblas
    BLASDEP := ../lib/librefblas.a
  endif
else
  $(error bad blas)
endif

# Compiler settings
DEBUGOPTS := -g
OPTIMOPTS := -O3
WARNOPTS  :=
AR := ar
ifeq ($(FC),gfortran)
  FOPTS       += -fimplicit-none -Wall
  WARN        := -std=f2003 -pedantic
  REFBLASWARN := -Wno-unused-dummy-argument -Wno-conversion
  TESTWARN    := -Wno-unused-dummy-argument -Wno-zerotrip
  DEBUGOPTS   += -fbacktrace -fcheck=all
  MOD         := -J
else ifeq ($(FC),pgfortran)#this one doesn't work with dispmodule
  MOD        := -module
else ifeq ($(FC),ifort)
  OPTIMOPTS  += -fast
  MOD        := -module
  AR         := xiar
else ifeq ($(FC),nagfor)
  FOPTS      += -u
  WARN       := -f2003 -info
  TESTWARN   := -f2003 -w=uda -info -w=unreffed -w=uparam
  REFBLASWARN:= -f2003 -w=uda -info -w=unreffed -w=obs
  DEBUGOPTS  += -gline -C=all
  OPTIMOPTS  := -O4
  MOD        := -mdir
else
  MOD        := -J #just guessing
endif

FOPTSREFBLAS := $(FOPTS) $(REFBLASWARN)
FOPTSTEST := $(FOPTS) $(TESTWARN)
FOPTS += $(WARN)

ifneq ($(debug),)
  FOPTS += $(DEBUGOPTS)
  LOPTS += $(DEBUGOPTS)
endif
ifneq ($(optimize),)
  FOPTS += $(OPTIMOPTS)
endif

FOPTSREFBLAS := $(FOPTS) $(REFBLASWARN)
FOPTS += $(WARNOPTS)
FOPTSTEST := $(FOPTS) $(TESTWARN)

# File name extensions:
OEXT       := .o
EEXT       :=

# Variables used by the non-recursive make:
ifeq ($(TOPLEVEL),true)
  LIB := lib
  SRC := src
  TOOLS := tools
else
  LIB := ../lib
  SRC := ../src
  TOOLS := ../tools
endif

endif #ifndef MAKE_INC
