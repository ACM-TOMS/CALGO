# Timing program for the Level 1 BLAS and selected LAPACK auxiliary
# routines.  This Makefile times the proposed new versions.
#
# E. Anderson
# May 9, 2016

ARCHFILE = ../etc/Makeinclude.${ARCH}
include $(ARCHFILE)
EXEDIR  = $(ARCH)
OBJDIR  = $(ARCH)
LIBDIR  = ../src/$(ARCH)
OUTDIR  = $(ARCH)

.SUFFIXES: .f .F .f90 .F90

BLASLIB = $(LIBDIR)/libblas1.a
SB1EXE = $(EXEDIR)/sb1tim
DB1EXE = $(EXEDIR)/db1tim
CB1EXE = $(EXEDIR)/cb1tim
ZB1EXE = $(EXEDIR)/zb1tim

AB1TEST = lsame.o lsamen.o iargcks.o la_xxvals.o
SB1TEST = slapy2.o slargv.o \
   stimmain.o sargcks.o sargckv.o sargckx.o srtg01.o second.o \
   stimamax.o stimasum.o stimaxpy.o stimcopy.o stimdot.o \
   stimnrm2.o stimpy2.o stimrot.o stimrtg.o stimscal.o stimswap.o \
   stomamax.o stomasum.o stomaxpy.o stomcopy.o stomdot.o \
   stomnrm2.o stompy2.o stomrot.o stomrtg.o stomscal.o stomswap.o
DB1TEST = dlapy2.o dlargv.o \
   dtimmain.o dargcks.o dargckv.o dargckx.o drtg01.o dsecnd.o \
   dtimamax.o dtimasum.o dtimaxpy.o dtimcopy.o dtimdot.o \
   dtimnrm2.o dtimpy2.o dtimrot.o dtimrtg.o dtimscal.o dtimswap.o \
   dtomamax.o dtomasum.o dtomaxpy.o dtomcopy.o dtomdot.o \
   dtomnrm2.o dtompy2.o dtomrot.o dtomrtg.o dtomscal.o dtomswap.o
CB1TEST = clargv.o \
   ctimmain.o sargcks.o cargcks.o cargckv.o cargckx.o crtg01.o second.o \
   ctimamax.o ctimasum.o ctimaxpy.o ctimcopy.o ctimdotc.o ctimdotu.o \
   ctimnrm2.o ctimrot.o ctimrtg.o ctimscal.o ctimsscal.o ctimsrot.o ctimswap.o \
   ctomamax.o ctomasum.o ctomaxpy.o ctomcopy.o ctomdotc.o ctomdotu.o \
   ctomnrm2.o ctomrot.o ctomrtg.o ctomscal.o ctomsscal.o ctomsrot.o ctomswap.o
ZB1TEST = zlargv.o \
   ztimmain.o dargcks.o zargcks.o zargckv.o zargckx.o zrtg01.o dsecnd.o \
   ztimamax.o ztimasum.o ztimaxpy.o ztimcopy.o ztimdotc.o ztimdotu.o \
   ztimnrm2.o ztimrot.o ztimrtg.o ztimscal.o ztimdscal.o ztimdrot.o ztimswap.o \
   ztomamax.o ztomasum.o ztomaxpy.o ztomcopy.o ztomdotc.o ztomdotu.o \
   ztomnrm2.o ztomrot.o ztomrtg.o ztomscal.o ztomdscal.o ztomdrot.o ztomswap.o

OBJECTS = $(addprefix $(OBJDIR)/,$(AB1TEST)) \
          $(addprefix $(OBJDIR)/,$(SB1TEST))
OBJECTD = $(addprefix $(OBJDIR)/,$(AB1TEST)) \
          $(addprefix $(OBJDIR)/,$(DB1TEST))
OBJECTC = $(addprefix $(OBJDIR)/,$(AB1TEST)) \
          $(addprefix $(OBJDIR)/,$(CB1TEST))
OBJECTZ = $(addprefix $(OBJDIR)/,$(AB1TEST)) \
          $(addprefix $(OBJDIR)/,$(ZB1TEST))

all: $(EXEDIR) $(OBJDIR) $(SB1EXE) $(DB1EXE) $(CB1EXE) $(ZB1EXE)

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

$(SB1EXE): $(OBJECTS)
	$(FC) -o $@ $(OBJECTS) $(MODOPT) $(OBJDIR) $(BLASLIB)

$(DB1EXE): $(OBJECTD)
	$(FC) -o $@ $(OBJECTD) $(MODOPT) $(OBJDIR) $(BLASLIB)

$(CB1EXE): $(OBJECTC)
	$(FC) -o $@ $(OBJECTC) $(MODOPT) $(OBJDIR) $(BLASLIB)

$(ZB1EXE): $(OBJECTZ)
	$(FC) -o $@ $(OBJECTZ) $(MODOPT) $(OBJDIR) $(BLASLIB)

check: $(SB1EXE) $(DB1EXE) $(CB1EXE) $(ZB1EXE)
	$(SB1EXE) < sbtime.in > $(OUTDIR)/sbtime.out
	$(DB1EXE) < dbtime.in > $(OUTDIR)/dbtime.out
	$(CB1EXE) < cbtime.in > $(OUTDIR)/cbtime.out
	$(ZB1EXE) < zbtime.in > $(OUTDIR)/zbtime.out

clean:
	/bin/rm -rf $(OBJECTS) $(OBJECTD) $(OBJECTC) $(OBJECTZ) \
        $(SB1EXE) $(DB1EXE) $(CB1EXE) $(ZB1EXE)

$(OBJDIR)/%.o: %.f $(ARCHFILE)
	$(FC) -c $(FFLAGS) -o $@ $(MODOPT) $(LIBDIR) $<

$(OBJDIR)/%.o: %.F $(ARCHFILE)
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $(MODOPT) $(LIBDIR) $<

$(OBJDIR)/%.o: %.f90 $(ARCHFILE)
	$(FC) -c $(FFLAGS) -o $@ $(MODOPT) $(LIBDIR) $<

$(OBJDIR)/%.o: %.F90 $(ARCHFILE)
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $(MODOPT) $(LIBDIR) $<
