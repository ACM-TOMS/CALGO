# Makefile for the Level 1 BLAS test suite
# Requires a file etc/Makeinclude.arch for this architecture
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
# April 25, 2016
#
ARCHFILE = ../etc/Makeinclude.${ARCH}
include $(ARCHFILE)
EXEDIR  = $(ARCH)
OBJDIR  = $(ARCH)
LIBDIR  = ../src/$(ARCH)
OUTDIR  = $(ARCH)

.SUFFIXES: .f .F .f90 .F90

BLASLIB = $(LIBDIR)/libblas1.a
SB1EXE = $(EXEDIR)/sb1tst
DB1EXE = $(EXEDIR)/db1tst
CB1EXE = $(EXEDIR)/cb1tst
ZB1EXE = $(EXEDIR)/zb1tst

ALLAUX = lsame.o lsamen.o iargcks.o
AB1TEST = la_xxvals.o
SB1TEST = sb1tmain.o sargcks.o sargckv.o sargckx.o \
   sb1tamax.o sb1tasum.o sb1taxpy.o sb1tcopy.o sb1tdot.o \
   sb1tnrm2.o sb1trot.o sb1trtg.o sb1tscal.o sb1tssq.o sb1tswap.o  \
   srtg01.o
DB1TEST = db1tmain.o dargcks.o dargckv.o dargckx.o \
   db1tamax.o db1tasum.o db1taxpy.o db1tcopy.o db1tdot.o \
   db1tnrm2.o db1trot.o db1trtg.o db1tscal.o db1tssq.o db1tswap.o  \
   drtg01.o
CB1TEST = cb1tmain.o cargcks.o cargckv.o cargckx.o sargcks.o \
   cb1tamax.o cb1tasum.o cb1taxpy.o cb1tcopy.o cb1tdotc.o cb1tdotu.o \
   cb1tnrm2.o cb1trot.o cb1trtg.o cb1tscal.o cb1tsrot.o cb1tsscal.o \
   cb1tssq.o cb1tswap.o crtg01.o slapy2.o
ZB1TEST = zb1tmain.o zargcks.o zargckv.o zargckx.o dargcks.o \
   zb1tamax.o zb1tasum.o zb1taxpy.o zb1tcopy.o zb1tdotc.o zb1tdotu.o \
   zb1tnrm2.o zb1trot.o zb1trtg.o zb1tscal.o zb1tdrot.o zb1tdscal.o \
   zb1tssq.o zb1tswap.o zrtg01.o dlapy2.o

AOBJS=$(addprefix $(OBJDIR)/,$(AB1TEST))
SOBJS=$(addprefix $(OBJDIR)/,$(SB1TEST)) \
      $(addprefix $(OBJDIR)/,$(ALLAUX))
DOBJS=$(addprefix $(OBJDIR)/,$(DB1TEST)) \
      $(addprefix $(OBJDIR)/,$(ALLAUX))
COBJS=$(addprefix $(OBJDIR)/,$(CB1TEST)) \
      $(addprefix $(OBJDIR)/,$(ALLAUX))
ZOBJS=$(addprefix $(OBJDIR)/,$(ZB1TEST)) \
      $(addprefix $(OBJDIR)/,$(ALLAUX))

all: $(ARCHFILE) $(EXEDIR) $(OBJDIR) $(SB1EXE) $(DB1EXE) $(CB1EXE) $(ZB1EXE)

# Exit if Makeinclude file does not exist
$(ARCHFILE):
	@echo "Check that the ARCH environment variable is set correctly"
	@exit 1

$(EXEDIR):
	mkdir -p $(EXEDIR)

ifneq ($(EXEDIR),$(OBJDIR))
$(OBJDIR):
	mkdir -p $(OBJDIR)
endif

$(SB1EXE): $(AOBJS) $(SOBJS) $(BLASLIB)
	$(FC) -o $@ $(AOBJS) $(SOBJS) $(BLASLIB)
$(DB1EXE): $(AOBJS) $(DOBJS) $(BLASLIB)
	$(FC) -o $@ $(AOBJS) $(DOBJS) $(BLASLIB)
$(CB1EXE): $(AOBJS) $(COBJS) $(BLASLIB)
	$(FC) -o $@ $(AOBJS) $(COBJS) $(BLASLIB)
$(ZB1EXE): $(AOBJS) $(ZOBJS) $(BLASLIB)
	$(FC) -o $@ $(AOBJS) $(ZOBJS) $(BLASLIB)

check: $(SB1EXE) $(DB1EXE) $(CB1EXE) $(ZB1EXE)
	$(SB1EXE) < sblas1t.in > $(OUTDIR)/sblas1t.out
	$(DB1EXE) < dblas1t.in > $(OUTDIR)/dblas1t.out
	$(CB1EXE) < cblas1t.in > $(OUTDIR)/cblas1t.out
	$(ZB1EXE) < zblas1t.in > $(OUTDIR)/zblas1t.out

clean:
	/bin/rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod \
           $(SB1EXE) $(DB1EXE) $(CB1EXE) $(ZB1EXE)

$(OBJDIR)/%.o:%.f $(ARCHFILE)
	$(FC) -c $(FFLAGS) -o $@ $(MODOPT) $(LIBDIR) $<

$(OBJDIR)/%.o:%.F $(ARCHFILE)
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $(MODOPT) $(LIBDIR) $<

$(OBJDIR)/%.o:%.f90 $(ARCHFILE)
	$(FC) -c $(FFLAGS) -o $@ $(MODOPT) $(LIBDIR) $<

$(OBJDIR)/%.o:%.F90 $(ARCHFILE)
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $(MODOPT) $(LIBDIR) $<
