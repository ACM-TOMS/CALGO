:: *** Windows DOS batch file to build a sample program using
:: ***              Talbot Suite DE
:: ***
:: *** author: Mariarosaria Rizzardi
:: ***         mariarosaria.rizzardi@uniparthenope.it
:: ***

@echo off

set TALBOTdir=../../TalbotSuiteDE
set EXdir=..

echo.
echo     *** Compile sources ***
@echo on
::gcc -c -std=gnu99 -Wall -pedantic -Wextra ...
:: with correction to NOPTS
::gcc -c -std=gnu99 -Wall -pedantic -Wno-unused-variable ./SEQ_main.c ./SEQ_LTsamples_fun.c %EXdir%/LTfun.c %EXdir%/ILTfun.c %EXdir%/LTsings.c %TALBOTdir%/COM/COM_Talbot_pack.c %TALBOTdir%/COM_DE/COM_Talbot_pack_DE.c %TALBOTdir%/FUN_DE/SEQ_Talbot_pack_DE.c %TALBOTdir%/FUN/SEQ_Talbot_pack.c
:: with no correction to NOPTS
gcc -c -std=gnu99 -Wall -pedantic -Wno-unused-variable ./SEQ_main.c ./SEQ_LTsamples_fun.c %EXdir%/LTfun.c %EXdir%/ILTfun.c %EXdir%/LTsings.c %TALBOTdir%/COM/COM_Talbot_pack.c %TALBOTdir%/COM_DE/COM_Talbot_pack_DE.c %TALBOTdir%/FUN_DE/SEQ_Talbot_pack_DE.c ./SEQ_Talbot_pack.c

@echo off
:: pause
echo.
echo.
echo.
echo     *** Build executable ***
echo.
@echo on
gcc -o ex0.exe .\*.o -lm

@echo off
:: pause
echo.
echo.
echo.
echo     *** Run executable as: ***
echo     ex0 n
echo.
echo     n = 1,2,3,4
echo     1)      20 t in [1000, 3000]
echo     2)     120 t in [1000, 3000]
echo     3)     120 t in [ 100,  500]
echo     4)      20 t in [ 100,  500]
echo.
echo     *** or execute: ***
echo     clean
echo     *** to remove object and executable files. ***
echo.
@echo on
