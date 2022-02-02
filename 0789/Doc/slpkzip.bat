rem Batchfile to create distribution archive of SLDRIVER package
rem It must live in the top level directory of the package
rem (usually C:\SL) and this must be the current directory when it
rem is called
rem It forms, in this same directory, a self-extracting PKZIP 2.05g
rem archive of the files in the package.
rem This is left in the file
rem      slzip.exe
rem in the directory referred to by %TEMP%.
rem If the lines at the end are uncommented it copies to drive A:
rem  - this zip file
rem  - the uncompressed READ.ME file
rem  - the PKZIP software (not neded for extraction but useful!)
rem !! Ensure disk in A: is blank to start with !!

rem !!! CONFIGURATION DATA !!!
set PK=C:\PK
set TEMP=C:\TEMP
set PATH=%PK%;%PATH%

rem Set 'archive' attribute on all files to be stored & on no others:
call slattrib

pkzip slzip -rp -i *.*
zip2exe slzip
del slzip.zip
move slzip.exe %TEMP%\slzip.exe

rem copy read.me a:
rem copy %TEMP%\slzip.exe a:\
rem xcopy %PK% a:\pk\ /s
