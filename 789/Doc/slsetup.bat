rem COPYRIGHT J D PRYCE 1998
rem PURPOSE
rem  Batchfile to do all the compilation & linking with Salford system
rem  To be adapted for other Fortran systems.
rem USAGE
rem  To be placed in the top-level directory of the package (usually
rem  c:\sl) and invoked from there.
rem NOTES
rem  This is provided as an alternative to the MAKEFILE but the latter
rem  is more convenient when it works.

rem ADJUST THE FOLLOWING CONFIGURATION SETTINGS AS NEEDED:
rem COPT = options for Salford FTN90 compiler
set COPT=/ca/debug

rem First we move to the main source-files directory. All compilation
rem is done there so that .MOD files are placed there.
cd sldriver
ftn90 ..\marcopak\marcomod.for%COPT%
cd ..\marcopak
link77 marcomod.lnk
cd ..\sldriver

ftn90 ..\sledge\sledgemd.for%COPT%
cd ..\sledge
link77 sledgemd.lnk
cd ..\sldriver

ftn90 ..\sleign\sleignmd.for%COPT%
cd ..\sleign
link77 sleignmd.lnk
cd ..\sldriver

ftn90 ..\sleign2\sleig2md.for%COPT%
cd ..\sleign2
link77 sleig2md.lnk
cd ..\sldriver

ftn90 slconsts.for%COPT%
ftn90 slpset.for%COPT%
ftn90 safeio.for%COPT%
ftn90 batchio.for%COPT%
ftn90 slutil.for%COPT%
ftn90 dbmod.for%COPT%
ftn90 sltstvar.for%COPT%
ftn90 ..\probsets\standard.for%COPT%
ftn90 ..\probsets\sample.for%COPT%
ftn90 sltstpak.for%COPT%
ftn90 solvrs.for%COPT%
ftn90 sldriver.for%COPT%
link77 standard.lnk
link77 sample.lnk
link77 bstandrd.lnk
link77 bsample.lnk
