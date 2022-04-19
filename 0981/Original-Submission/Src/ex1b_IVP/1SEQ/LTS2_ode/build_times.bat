cls
@echo off

set TALBOTdir=../../../TalbotSuiteDE
set EXdir=../..

echo.
echo BUILD EXECUTABLE FOR TIME TEST
echo     [build_times.bat]
echo.
echo.
echo Example folder is: "%EXdir%"
echo Talbot Suite folder is: "%TALBOTdir%"

echo.
echo     *** Compile sources ***
@echo on
gcc -c -std=gnu99 -Wall -pedantic -Wextra ./SEQ_main_TIMES.c ./SEQ_LTsamples_ode.c ./ode.c %EXdir%/LTsings2.c %TALBOTdir%/COM/COM_Talbot_pack.c %TALBOTdir%/COM_DE/COM_Talbot_pack_DE.c %TALBOTdir%/FUN_DE/SEQ_Talbot_pack_DE.c

@echo off
:: pause
echo.
echo.
echo.
echo     *** Build executable ***
echo.
@echo on
gcc -o ex1_time.exe ./SEQ_main_TIMES.o ./SEQ_LTsamples_ode.o ./ode.o ./LTsings2.o ./COM_Talbot_pack.o ./COM_Talbot_pack_DE.o ./SEQ_Talbot_pack_DE.o -lm

@echo off
:: pause
echo.
echo.
echo.
echo     *** Run executable as: ***
echo     ex1_time tol
echo.
echo     where:
echo.
echo     tol = 1e-6, 1e-8, 1e-10, 1e-12
echo.
echo.
echo     *** or execute: ***
echo     clean
echo     *** to remove object and executable files. ***
echo.
@echo on

:: OUTPUT TO THE DISPLAY
::ex1_time 1e-6

:: OUTPUT TO A FILE
::@echo off
::echo.
::echo     *** Output to file: out.txt ***
::echo.
::@echo on
::ex1_time 1e-6 > out.txt
