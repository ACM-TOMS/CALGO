:: *** Windows DOS batch file to build and run a sample program using
:: ***              Talbot Suite DE
:: ***
:: *** author: Mariarosaria Rizzardi
:: ***         mariarosaria.rizzardi@uniparthenope.it
:: ***

@echo off
cls

set TALBOTdir=../../TalbotSuiteDE
set EXdir=..

echo Example folder is: "%EXdir%"
echo Talbot Suite folder is: "%TALBOTdir%"

echo.
echo     *** Compile sources ***
@echo on
::gcc -c -std=gnu99 -Wall -pedantic -Wextra ./SEQ_main.c ./SEQ_LTsamples_fun.c %EXdir%/LTfun.c %EXdir%/ILTfun.c %EXdir%/LTsings.c %TALBOTdir%/COM/COM_Talbot_pack.c %TALBOTdir%/COM_DE/COM_Talbot_pack_DE.c %TALBOTdir%/FUN/SEQ_Talbot_pack.c %TALBOTdir%/FUN_DE/SEQ_Talbot_pack_DE.c
gcc -c -std=gnu99 -Wall -pedantic -Wno-unused-variable ./SEQ_main.c ./SEQ_LTsamples_fun.c %EXdir%/LTfun.c %EXdir%/ILTfun.c %EXdir%/LTsings.c %TALBOTdir%/COM/COM_Talbot_pack.c %TALBOTdir%/COM_DE/COM_Talbot_pack_DE.c %TALBOTdir%/FUN/SEQ_Talbot_pack.c %TALBOTdir%/FUN_DE/SEQ_Talbot_pack_DE.c

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
echo.
echo ********************************************
echo.
echo    Choose a test:
echo    1)      20 t in [1000, 3000]
echo    2)     120 t in [1000, 3000]
echo    3)     120 t in [ 100,  500]
echo    4)      20 t in [ 100,  500]
echo.
echo ********************************************
echo.
set /p choice=Enter [1..4]?: 
echo.
echo choice = %choice%
echo.
echo.
if %choice% LSS 1 (
    echo ********************************************
    echo.
    echo     *** Error: entered a value less than 1 ***
    echo         EXIT!
    echo.
    echo ********************************************
    echo.
    pause
    exit
)
if %choice% GTR 4 (
    echo ********************************************
    echo.
    echo     *** Error: entered a value greater than 4 ***
    echo         EXIT!
    echo.
    echo ********************************************
    echo.
    pause
    exit
)

:: choice is between 1 and 4
echo ********************************************
echo.
echo     *** Run executable as: ***
echo         ex0.exe %choice%
echo.
echo ********************************************
echo.
@echo on
:: OUTPUT TO THE DISPLAY (1,2 or 3)
ex0.exe %choice%


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

:: OUTPUT TO A FILE
::@echo off
::echo.
::echo     *** Output to file: out.txt ***
::echo.
::@echo on
::ex0.exe 1 > out.txt
