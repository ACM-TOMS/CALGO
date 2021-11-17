# File:Makefile, revision of Aug 1998
# Purpose:
# Manage the files making up the SLDRIVER package, when using the
# NAG f90 system under Unix
# Author : John Pryce, RMCS, Shrivenham, Swindon, UK
#(pryce@rmcs.cranfield.ac.uk)
# Disclaimer:
# This file is provided as is, with no claim that it operates
# correctly in all circumstances.
# References:
# SLDRIVER User Guide & Tutorial, RMCS Tech. Rep. SEAS/CISE/96/JDP01
#
###################### INSTRUCTIONS FOR USE ##############################
#1.This file must be placed in, and invoked from, the 'sldriver' directory
#under the home directory of the SL package. This ensures that the
#'.mod' (module specification) files are put in this directory.
#
#2.Set SL to the home directory of the SL package:
SL=/nfs/quince/d2/home/cur/trh/SandPit/CALGO/789
#
#3.Define the compiler
F90=f90 -c
LINK=f90
#4.If desired, reset these compilation options:
OPTS=
#5.The extension of a module specification file:
MOD=mod
#The following can be typed at the DOS prompt:
#COMMAND |EFFECT
#------------------------------------------------------------------------
#mk|makes standard.x if necessary, and launches it.
#mk standard.x |as above but without launching it.
#mk sample |makes sample.x if necessary, and launches it.
#mk bstandrd.x |makes bstandrd.x if necessary, doesn't launch it.
#mk sample.x |as above but without launching it.
#etc... (the above are the most useful).
#########################################################################

## Run the main executable, which uses the STANDARD problem set
standard: standard.x
	standard.x
standard.x: sldriver.o slutil.o solvrs.o dbmod.o \
	safeio.o sltstpak.o standard.o
	$(LINK) -o standard.x sldriver.o slutil.o solvrs.o dbmod.o safeio.o \
	sltstvar.o slpset.o \
	sledgemd.o sleignmd.o sleig2md.o marcomod.o \
	sltstpak.o standard.o

## Alternative executable, using SAMPLE problem set
sample: sample.x
	sample.x
sample.x: sldriver.o slutil.o solvrs.o dbmod.o \
	safeio.o sltstpak.o sample.o
	$(LINK) -o sample.x sldriver.o slutil.o solvrs.o dbmod.o safeio.o \
	sltstvar.o slpset.o \
	sledgemd.o sleignmd.o sleig2md.o marcomod.o \
	sltstpak.o sample.o

## Batch-run executable, using STANDARD problem set
bstandrd.x: sldriver.o slutil.o solvrs.o dbmod.o \
	batchio.o sltstpak.o standard.o
	$(LINK) -o bstandrd.x sldriver.o slutil.o solvrs.o dbmod.o batchio.o \
	sltstvar.o slpset.o \
	sledgemd.o sleignmd.o sleig2md.o marcomod.o \
	sltstpak.o standard.o

## Batch-run executable, using SAMPLE problem set
bsample.x: sldriver.o slutil.o solvrs.o dbmod.o \
	batchio.o sltstpak.o sample.o
	$(LINK) -o bsample.x sldriver.o slutil.o solvrs.o dbmod.o batchio.o \
	sltstvar.o slpset.o \
	sledgemd.o sleignmd.o sleig2md.o marcomod.o \
	sltstpak.o sample.o

##Make all the executables:
all.x:  standard.x bstandrd.x sample.x bsample.x

## General compilation rule
.f.o:
	$(F90) $(OPTS) $<

## The main-program source
sldriver.o: sldriver.f dbmod.o safeio.o slutil.o solvrs.o \
	slconsts.o slpset.o sltstpak.o
	$(F90) $(OPTS) sldriver.f

## The database facilities
dbmod.o:dbmod.f safeio.o slutil.o slconsts.o slpset.o
	$(F90) $(OPTS) dbmod.f

## The solver interfaces
solvrs.o: solvrs.f safeio.o slutil.o slconsts.o \
	sltstpak.o sledge.o sleign.o marcopak.o sleign2.o
	$(F90) $(OPTS) solvrs.f

## Utilities
#safeio.o: safeio.f
#$(F90) $(OPTS) safeio.f
#batchio.o:batchio.f
#$(F90) $(OPTS) batchio.f
#slutil.o: slutil.f
#$(F90) $(OPTS) slutil.f
#slconsts.o: slconsts.f
#$(F90) $(OPTS) slconsts.f
#slpset.o: slpset.f
#$(F90) $(OPTS) slpset.f

## SLTSTPAK: problem/solver interface routines
sltstpak.o: sltstpak.f sltstvar.o testmod.$(MOD)
	$(F90) $(OPTS) sltstpak.f

## Problem-set routines
##BUG HERE! I don't know how to make MK distinguish between different
##problem sets, so at present it only seems to work if one makes the
##standard executable before the sample one (TESTMOD.MOD gets out of
##kilter otherwise?)
testmod.mod: standard.o
##testmod.mod: sample.o

standard.o: $(SL)/probsets/standard.f sltstvar.o
	$(F90) $(OPTS) $(SL)/probsets/standard.f
sample.o:   $(SL)/probsets/sample.f sltstvar.o
	$(F90) $(OPTS) $(SL)/probsets/sample.f

## The solvers. They live in their own directories.
## Compilation is to be done from within the SLDRIVER directory
## so that the .o and .mod files get stored there.
sledge.o: $(SL)/sledge/sledgemd.f
	$(F90) $(OPTS) $(SL)/sledge/sledgemd.f

sleign.o: $(SL)/sleign/sleignmd.f
	$(F90) $(OPTS) $(SL)/sleign/sleignmd.f

marcopak.o: $(SL)/marcopak/marcomod.f
	$(F90) $(OPTS) $(SL)/marcopak/marcomod.f

sleign2.o:$(SL)/sleign2/sleig2md.f
	$(F90) $(OPTS) $(SL)/sleign2/sleig2md.f
