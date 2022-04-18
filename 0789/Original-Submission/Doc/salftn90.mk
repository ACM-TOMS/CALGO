# File: Makefile for SLDRIVER package, revision of Aug 1998
#!Using LINK77 which needs .lnk files (instead of SLINK)
# Purpose:
# Manage the files making up the SLDRIVER package, when using the
# Salford FTN90 system under DOS/Windows 3.1 or Windows 95
# Author : John Pryce, RMCS, Shrivenham, Swindon, UK
#(pryce@rmcs.cranfield.ac.uk)
# Disclaimer:
# This file is provided as is, with no claim that it operates
# correctly in all circumstances.
# References:
# SLDRIVER User Guide & Tutorial, RMCS Tech. Rep. SEAS/CISE/96/JDP01
#
###################### INSTRUCTIONS FOR USE ##############################
# This file must be placed in, and invoked from, the 'sldriver' directory
#under the home directory of the SL package. This ensures that the
#'.$M' (module specification) files are put in this directory.
#2.Set SL to the home directory of the SL package:
SL=C:\sl4.1
BIN=$(SL)\SLDRIVER
#
#3.Define the compiler & linking mechanism
F90=FTN90
LINK=LINK77
#4.If desired, reset these compilation options.
#  /BINARY stuff is needed so object file goes in current directory
#  even if source isn't.
#OPTS=/CHECK /SILENT /ERRORLOG /BINARY $(BIN)\$@
OPTS=/CHECK_ALL /SILENT /BINARY $(BIN)\$@
#5.Various file extensions:
F=for # Fortran (fixed format) file
M=mod # module specification file:
O=obj # object file
X=exe # executable
.SUFFIXES: .$F .$O
## General compilation rule
.$F.$O:
	$(F90) $(OPTS) $<

#The following can be typed at the DOS prompt:
#COMMAND |EFFECT
#------------------------------------------------------------------------
#mk              |makes standard.exe if necessary, and launches it.
#mk standard.exe |as above but without launching it.
#mk sample       |makes sample.exe if necessary, and launches it.
#mk bstandrd.exe |makes bstandrd.exe if necessary, doesn't launch it.
#mk sample.exe   |as above but without launching it.
#etc... (the above are the most useful).
#########################################################################

## Run the main executable, which uses the STANDARD problem set
standard: standard.$X
        RUN77 standard.$X
standard.$X: sldriver.$O slutil.$O solvrs.$O dbmod.$O \
        safeio.$O sltstpak.$O standard.$O
#        $(LINK) -out:standard.$X sldriver.$O slutil.$O solvrs.$O dbmod.$O safeio.$O \
#        sltstvar.$O slpset.$O \
#        sledge.$O sleign.$O sleign2.$O marcopak.$O \
#        sltstpak.$O standard.$O
         LINK77 standard.lnk
## Alternative executable, using SAMPLE problem set
sample: sample.$X
        RUN77 sample.$X
sample.$X: sldriver.$O slutil.$O solvrs.$O dbmod.$O \
        safeio.$O sltstpak.$O sample.$O
#        $(LINK) -out:sample.$X sldriver.$O slutil.$O solvrs.$O dbmod.$O safeio.$O \
#        sltstvar.$O slpset.$O \
#        sledge.$O sleign.$O sleign2.$O marcopak.$O \
#        sltstpak.$O sample.$O
         LINK77 sample.lnk

## Batch-run executable, using STANDARD problem set
bstandrd.$X: sldriver.$O slutil.$O solvrs.$O dbmod.$O \
        batchio.$O sltstpak.$O standard.$O
#        $(LINK) -out:standard.$X sldriver.$O slutil.$O solvrs.$O dbmod.$O batchio.$O \
#        sltstvar.$O slpset.$O \
#        sledge.$O sleign.$O sleign2.$O marcopak.$O \
#        sltstpak.$O standard.$O
         LINK77 bstandrd.lnk

## Batch-run executable, using SAMPLE problem set
bsample.$X: sldriver.$O slutil.$O solvrs.$O dbmod.$O \
        batchio.$O sltstpak.$O sample.$O
#        $(LINK) -out:sample.$X sldriver.$O slutil.$O solvrs.$O dbmod.$O batchio.$O \
#        sltstvar.$O slpset.$O \
#        sledge.$O sleign.$O sleign2.$O marcopak.$O \
#        sltstpak.$O sample.$O
         LINK77 bsample.lnk

## The main-program source
sldriver.$O: sldriver.$F dbmod.$O safeio.$O slutil.$O solvrs.$O \
        slconsts.$O slpset.$O sltstpak.$O
        $(F90) $(OPTS) sldriver.$F

## The database facilities
dbmod.$O:dbmod.$F safeio.$O slutil.$O slconsts.$O slpset.$O
        $(F90) $(OPTS) dbmod.$F

## The solver interfaces
solvrs.$O: solvrs.$F safeio.$O slutil.$O slconsts.$O \
        sltstpak.$O sledge.$O sleign.$O marcopak.$O sleign2.$O
        $(F90) $(OPTS) solvrs.$F

## Utilities
#safeio.$O: safeio.$F
#$(F90) $(OPTS) safeio.$F
#batchio.$O:batchio.$F
#$(F90) $(OPTS) batchio.$F
#slutil.$O: slutil.$F
#$(F90) $(OPTS) slutil.$F
#slconsts.$O: slconsts.$F
#$(F90) $(OPTS) slconsts.$F
#slpset.$O: slpset.$F
#$(F90) $(OPTS) slpset.$F

## SLTSTPAK: problem/solver interface routines
sltstpak.$O: sltstpak.$F sltstvar.$O testmod.$M
        $(F90) $(OPTS) sltstpak.$F

## Problem-set routines
##BUG HERE! I don't know how to make MK distinguish between different
##problem sets, so at present it only seems to work if one changes this
##rule & deletes standard.$O & sample.$O before running 'make' to
##create the executable for a different problem set.
##(TESTMOD.MOD gets out of kilter otherwise?)
testmod.$M: standard.$O
##testmod.$M: sample.$O

standard.$O: $(SL)\PROBSETS\standard.$F sltstvar.$O
        $(F90) $(OPTS) $(SL)\PROBSETS\standard.$F
sample.$O:   $(SL)\PROBSETS\sample.$F sltstvar.$O
        $(F90) $(OPTS) $(SL)\PROBSETS\sample.$F

## The solvers. They live in their own directories.
## Compilation is to be done from within the SLDRIVER directory
## so that the .$O and .$M files get stored there.
sledge.$O: $(SL)\SLEDGE\sledgemd.$F
        $(F90) $(OPTS) $(SL)\SLEDGE\sledgemd.$F

sleign.$O: $(SL)\SLEIGN\sleignmd.$F
        $(F90) $(OPTS) $(SL)\SLEIGN\sleignmd.$F

marcopak.$O: $(SL)\MARCOPAK\marcomod.$F
        $(F90) $(OPTS) $(SL)\MARCOPAK\marcomod.$F

sleign2.$O:$(SL)\SLEIGN2\sleig2md.$F
        $(F90) $(OPTS) $(SL)\SLEIGN2\sleig2md.$F
