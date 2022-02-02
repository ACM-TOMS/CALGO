cls
@echo off

set TALBOTdir=../../../TalbotSuiteDE
set EXdir=../..

echo Example folder is: "%EXdir%"
echo Talbot Suite folder is: "%TALBOTdir%"

echo.
echo     *** Compile sources ***
@echo on

gfortran -c twpbvp.f prtfun.f testfun.f

gcc -c -std=gnu99 -Wall -pedantic -Wextra ./SEQ_main_ACCURACY.c ./SEQ_LTsamples_twpbvp.c %EXdir%/ILTfun2.c %EXdir%/LTsings2.c %TALBOTdir%/COM/COM_Talbot_pack.c %TALBOTdir%/COM_DE/COM_Talbot_pack_DE.c %TALBOTdir%/FUN_DE/SEQ_Talbot_pack_DE.c

@echo off
:: pause
echo.
echo.
echo.
echo     *** Build executable ***
echo.
@echo on
gcc -o ex3_acc.exe ./SEQ_main_ACCURACY.o ./SEQ_LTsamples_twpbvp.o ./ILTfun2.o ./LTsings2.o ./COM_Talbot_pack.o ./COM_Talbot_pack_DE.o ./SEQ_Talbot_pack_DE.o twpbvp.o prtfun.o testfun.o -lm -lgfortran

@echo off
:: pause
echo.
echo.
echo.
echo     *** Run executable as: ***
echo     ex3_acc tol
echo.
echo     where:
echo.
echo     tol = 1e-6, 1e-8, 1e-10, 1e-12
echo.
@echo on

:: OUTPUT TO THE DISPLAY
::ex3_acc 1e-6

:: OUTPUT TO A FILE
::@echo off
::echo.
::echo     *** Output to file: out.txt ***
::echo.
::@echo on
::ex3_acc 1e-6 > out.txt
