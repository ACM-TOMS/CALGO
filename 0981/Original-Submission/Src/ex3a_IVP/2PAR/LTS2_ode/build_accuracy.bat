cls
@echo off
:: *** Windows DOS batch file to build and run a sample program using
:: ***              Talbot Suite DE
:: ***
:: *** author: Mariarosaria Rizzardi
:: ***         mariarosaria.rizzardi@uniparthenope.it
:: ***
:: *** To run:	build_accuracy

set TALBOTdir=../../../TalbotSuiteDE
set EXdir=../..

echo Example folder is: "%Exdir%"
echo Talbot Suite folder is: "%TALBOTdir%"

echo.
echo     *** Compile sources ***
@echo on
::gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wextra ./OMP_main_ACCURACY.c ./OMP_LTsamples_ode.c ./ode.c %EXdir%/ILTfun2.c %EXdir%/LTsings2.c %TALBOTdir%/COM/COM_Talbot_pack.c %TALBOTdir%/COM_DE/COM_Talbot_pack_DE.c %TALBOTdir%/FUN_DE/OMP_Talbot_pack_DE.c
gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wno-unused-but-set-variable -Wno-comment ./OMP_main_ACCURACY.c ./OMP_LTsamples_ode.c ./ode.c %EXdir%/ILTfun2.c %EXdir%/LTsings2.c %TALBOTdir%/COM/COM_Talbot_pack.c %TALBOTdir%/COM_DE/COM_Talbot_pack_DE.c %TALBOTdir%/FUN_DE/OMP_Talbot_pack_DE.c

@echo off
:: pause
echo.
echo.
echo.
echo     *** Build executable ***
echo.
@echo on
gcc -o ex3_acc.exe ./OMP_main_ACCURACY.o ./OMP_LTsamples_ode.o ./ode.o ./ILTfun2.o ./LTsings2.o ./COM_Talbot_pack.o ./COM_Talbot_pack_DE.o ./OMP_Talbot_pack_DE.o -lm -lgomp

@echo off
:: pause
echo.
echo.
echo.
echo     *** Run executable as: ***
echo.
echo     ex3_acc tol thrds1
echo     or
echo     ex3_acc tol thrds1 thrds2
echo.
echo     where:
echo.
echo     tol = tolerance (1e-6, 1e-8, 1e-10, 1e-12)
echo.
echo     thrds1, thrds2 number of Open MP threads
echo                    (both for nested parallelism such that
echo                    thrds1*thrds2 is the total number of threads)
echo     If only thrds1 is given, then thrds2 is set to 1.
echo.
@echo on
