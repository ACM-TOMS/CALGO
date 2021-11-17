cls
@echo off

set TALBOTdir=../../../TalbotSuiteDE
set EXdir=../..

echo Example folder is: "%Exdir%"
echo Talbot Suite folder is: "%TALBOTdir%"

echo.
echo     *** Compile sources ***
@echo on

::the extension has been changed into .F to abilitate the preprocessing directives
gfortran -c -fopenmp twpbvp.F prtfun.F testfun.F

::gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wextra ./OMP_main_TIMES.c ./OMP_LTsamples_twpbvp.c %EXdir%/LTsings2.c %TALBOTdir%/COM/COM_Talbot_pack.c %TALBOTdir%/COM_DE/COM_Talbot_pack_DE.c %TALBOTdir%/FUN_DE/OMP_Talbot_pack_DE.c
gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wno-unused-but-set-variable -Wno-comment ./OMP_main_TIMES.c ./OMP_LTsamples_twpbvp.c %EXdir%/LTsings2.c %TALBOTdir%/COM/COM_Talbot_pack.c %TALBOTdir%/COM_DE/COM_Talbot_pack_DE.c %TALBOTdir%/FUN_DE/OMP_Talbot_pack_DE.c

@echo off
:: pause
echo.
echo.
echo.
echo     *** Build executable ***
echo.
@echo on
gcc -o ex4_time.exe ./OMP_main_TIMES.o ./OMP_LTsamples_twpbvp.o ./LTsings2.o ./COM_Talbot_pack.o ./COM_Talbot_pack_DE.o ./OMP_Talbot_pack_DE.o twpbvp.o prtfun.o testfun.o -lm -lgfortran -lgomp

@echo off
:: pause
echo.
echo.
echo.
echo     *** Run executable as: ***
echo.
echo     ex4_time tol jFUN NTval NXval
echo.
echo     where:
echo.
echo     tol  = tolerance (1e-6, 1e-8, 1e-10, 1e-12)
echo.
echo     jFUN = 1 (coarse-grain parallelism)
echo          = 2 ( fine-grain  parallelism)
echo          = 3 (   nested    parallelism)
echo.
echo     NTval: number of t-values in [Tmin, Tmax]
echo.
echo     NXval: number of x-values in [Xmin, Xmax]
echo.
@echo on

:: OUTPUT TO A FILE
::@echo off
::echo.
::echo     *** Output to file: out.txt ***
::echo.
::@echo on
::ex4_time 1e-6 1 120 5 > out.txt
