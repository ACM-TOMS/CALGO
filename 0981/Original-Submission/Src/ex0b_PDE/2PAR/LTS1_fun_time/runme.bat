cls
@echo off

set TALBOTdir=../../../TalbotSuiteDE
set EXdir=../..

echo.
echo BUILD EXECUTABLE FOR TIME TEST
echo       Example 0b [Duffy] - LT samples by a function
echo              Test on execution times
echo Comparison between Talbot Suite and Talbot Suite DE
echo              OpenMP parallel versions
echo.
echo.
echo Example folder is: "%EXdir%"
echo Talbot Suite folder is: "%TALBOTdir%"

echo.
echo     *** Compile sources ***
@echo on
::gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wextra -Wno-unused-variable ...
gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wno-unused-but-set-variable -Wno-unused-parameter ./OMP_main_TIMES.c ./OMP_LTsamples_fun.c %EXdir%/LTfun2.c %EXdir%/LTsings2.c %TALBOTdir%/COM/COM_Talbot_pack.c %TALBOTdir%/COM_DE/COM_Talbot_pack_DE.c %TALBOTdir%/FUN_DE/OMP_Talbot_pack_DE.c

@echo off
:: pause
echo.
echo.
echo.
echo     *** Build executable ***
echo.
@echo on
gcc -o ex0.exe ./OMP_main_TIMES.o ./OMP_LTsamples_fun.o ./LTfun2.o ./LTsings2.o ./COM_Talbot_pack.o ./COM_Talbot_pack_DE.o ./OMP_Talbot_pack_DE.o -lm -lgomp

@echo off
:: pause
echo.
echo.
echo.
echo     *** Run executable as: ***
echo     ex0 1e-12 1 120 120
echo.
@echo on
ex0 1e-12 1 120 120

@echo off
:: pause
echo.
echo.
echo.
echo     *** Run executable as: ***
echo     ex0 tol jFUN NTval NXval
echo.
echo     where:
echo           tol  :  tolerance (1e-6, 1e-8, 1e-10, 1e-12)
echo           jFUN :  number of OMP function to be tested
echo                  = 1 OMP_Talbot11 [coarse-grain parallelism for modified method]
echo                  = 2 OMP_Talbot12 [ fine-grain  parallelism for modified method]
echo                  = 3 OMP_Talbot13 [   nested    parallelism for modified method]
echo           NTval:  number of t-values {5, 20, 120}
echo           NXval:  number of x-values {5, 20, 120}
echo.
echo.
echo.
echo     *** or execute: ***
echo     clean
echo     *** to remove object and executable files. ***
echo.
@echo on

:: OUTPUT TO THE DISPLAY
::ex0 1e-6 1 20 20

:: OUTPUT TO A FILE
::@echo off
::echo.
::echo     *** Output to file: out.txt ***
::echo.
::@echo on
::ex0 1e-6 1 20 20 > out.txt
