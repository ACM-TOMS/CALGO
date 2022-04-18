rem Batchfile to create Unix format distrib archive of SLDRIVER package
rem It forms a new directory structure under \SLUNIX holding the
rem DOS2UNIX'ed versions of the files, which it then TAR's and GZIP's
rem into a file SLTAR.Z in the main SL directory.
rem Uses the DOS versions of Unix tools TAR GZIP etc., must be in Path!

rem *********CHANGE THIS CONFIGURATION INFO IF NECESSARY***************
set SU=C:\slunix
set SL=C:\sl4.1

rem Set 'archive' attribute on all files to be stored & on no others:
cd %SL%
call slattrib

deltree %SU%
md %SU%
rem form copy of all these files
xcopy %SL% %SU%\ /a/s
cd %SU%
for %%f in (*.*) do call d2u %%f
cd %SU%\d02kef
for %%f in (*.*) do call d2u %%f
cd %SU%\marcopak
for %%f in (*.*) do call d2u %%f
cd %SU%\sledge
for %%f in (*.*) do call d2u %%f
cd %SU%\sleign
for %%f in (*.*) do call d2u %%f
cd %SU%\sleign2
for %%f in (*.*) do call d2u %%f
cd %SU%\sleign2\extras
for %%f in (*.*) do call d2u %%f
cd %SU%\sldriver
for %%f in (*.*) do call d2u %%f
cd %SU%\sldriver\standard
for %%f in (*.*) do call d2u %%f
cd %SU%\sldriver\standard\truevals
for %%f in (*.*) do call d2u %%f
cd %SU%\sldriver\sample
for %%f in (*.*) do call d2u %%f
cd %SU%\sldriver\sample\truevals
for %%f in (*.*) do call d2u %%f
cd %SU%\probsets
for %%f in (*.*) do call d2u %%f
cd %SU%\latex
for %%f in (*.*) do call d2u %%f

rem all dos2unix'ed, now archive them
cd %SU%
tar -f %SL%\sl.tar -cv .
cd %SL%
gzip -c sl.tar > sltar.z
