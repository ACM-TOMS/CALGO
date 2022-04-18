:: *** Windows DOS batch file to build a sample program using
:: ***              Talbot Suite DE
:: ***
:: *** author: Mariarosaria Rizzardi
:: ***         mariarosaria.rizzardi@uniparthenope.it
:: ***

cls
@echo off

set TALBOTdir=../../TalbotSuiteDE
set EXdir=..

echo Example folder is: "%EXdir%"
echo Talbot Suite folder is: "%TALBOTdir%"

echo.
echo     *** Compile sources ***
@echo on
::gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wextra ...
gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wno-unused-variable -Wno-comment ./OMP_main.c ./SEQ_LTsamples_fun.c ./OMP_LTsamples_fun.c ./LTfun2.c %EXdir%/LTfun.c %EXdir%/ILTfun.c %EXdir%/LTsings.c %TALBOTdir%/COM/COM_Talbot_pack.c %TALBOTdir%/COM_DE/COM_Talbot_pack_DE.c %TALBOTdir%/FUN_DE/OMP_Talbot_pack_DE.c %TALBOTdir%/FUN_DE/SEQ_Talbot_pack_DE.c

@echo off
:: pause
echo.
echo.
echo.
echo     *** Build executable ***
echo.
@echo on
gcc -o ex0.exe ./*.o -lm -lgomp


@echo off
echo.
echo.
echo.
echo     *** Run executable as: ***
echo     ex0 3 1
echo.
@echo on
ex0 3 1


@echo off
echo.
echo.
echo.
echo     *** Run executable as: ***
echo     ex0 n jFUN
echo.
echo.    where
echo        n = 1,2 or 3
echo        1)      20 t in [1000, 3000]    (both times and errors if PRINTflag=true)
echo        2)     120 t in [1000, 3000]    (only times)
echo        3)     120 t in [ 100,  500]    (only times)
echo.
echo        jFUN = 0 : SEQ_Talbot1_DE  - sequential modified Talbot's method
echo             = 1 : OMP_Talbot11_DE - OpenMP coarse-grained parallel modified Talbot's method
echo             = 2 : OMP_Talbot12_DE - OpenMP fine-grained parallel modified Talbot's method
echo             = 3 : OMP_Talbot13_DE - OpenMP nested coarse/fine-grained parallel modified Talbot's method
echo             = 4 : SEQ_Talbot2_DE  - sequential classical Talbot's method
echo.
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
::ex0.exe 3 1 > out.txt
