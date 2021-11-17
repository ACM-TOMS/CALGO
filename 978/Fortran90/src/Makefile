# Makefile for the Level 1 BLAS library
# Requires a file ../etc/Makeinclude.arch for this architecture
# Some sample Makeinclude files are provided in the etc directory.
#
# The environment variable ARCH must be set to match, for example
#    setenv ARCH intel
# The object files, module files, executable files, and output
# files are all created in the ARCH directory by default.
#
# This organization allows you to test multiple architectures
# without overwriting the files created in previous tests.
#
# Edward Anderson
# May 5, 2016
#
ARCHFILE = ../etc/Makeinclude.${ARCH}
include $(ARCHFILE)
LIBDIR  = $(ARCH)
OBJDIR  = $(ARCH)

.SUFFIXES: .f .F .f90 .F90

BLASLIB = $(LIBDIR)/libblas1.a

FMODS = la_constants32.o
FMODD = la_constants.o
FMOD  = $(FMODS) $(FMODD) la_xisnan.o
SBLAS1  = isamax.o \
   sasum.o saxpy.o scopy.o sdot.o slartg.o slassq.o snrm2.o \
   srot.o srotg.o sscal.o sswap.o
DBLAS1  = idamax.o \
   dasum.o daxpy.o dcopy.o ddot.o dlartg.o dlassq.o dnrm2.o \
   drot.o drotg.o dscal.o dswap.o
CBLAS1  = icamax.o \
   scasum.o caxpy.o ccopy.o cdotc.o cdotu.o clartg.o classq.o scnrm2.o \
   crot.o csrot.o crotg.o csscal.o cscal.o cswap.o
ZBLAS1  = izamax.o \
   dzasum.o zaxpy.o zcopy.o zdotc.o zdotu.o zlartg.o zlassq.o dznrm2.o \
   zrot.o zdrot.o zrotg.o zdscal.o zscal.o zswap.o

BLOBJS=$(addprefix $(OBJDIR)/,$(FMOD)) \
         $(addprefix $(OBJDIR)/,$(SBLAS1)) \
         $(addprefix $(OBJDIR)/,$(DBLAS1)) \
         $(addprefix $(OBJDIR)/,$(CBLAS1)) \
         $(addprefix $(OBJDIR)/,$(ZBLAS1))

all: $(LIBDIR) $(OBJDIR) $(BLASLIB)

$(LIBDIR):
	mkdir -p $(LIBDIR)

ifneq ($(LIBDIR),$(OBJDIR))
$(OBJDIR):
	mkdir -p $(OBJDIR)
endif

$(BLASLIB): $(BLOBJS)
	ar cr $(BLASLIB) $(BLOBJS)

clean:
	/bin/rm -rf $(BLASLIB) $(OBJDIR)/*.o $(OBJDIR)/*.mod

$(OBJDIR)/%.o: %.f $(ARCHFILE)
	$(FC) -c $(FFLAGS) -o $@ $(MODOPT) $(OBJDIR) $<

$(OBJDIR)/%.o: %.F $(ARCHFILE)
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $(MODOPT) $(OBJDIR) $<

$(OBJDIR)/%.o: %.f90 $(ARCHFILE)
	$(FC) -c $(FFLAGS) -o $@ $(MODOPT) $(OBJDIR) $<

$(OBJDIR)/%.o: %.F90 $(ARCHFILE)
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $(MODOPT) $(OBJDIR) $<
