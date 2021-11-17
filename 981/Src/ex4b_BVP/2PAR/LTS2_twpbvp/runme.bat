:: *** Windows DOS batch file to build and run a sample program using
:: ***              Talbot Suite DE
:: ***
:: *** author: Mariarosaria Rizzardi
:: ***         mariarosaria.rizzardi@uniparthenope.it
:: ***

cls
@echo off
echo.
echo    Select:
echo    1)  accuracy test
echo    2)  time test 
echo.
set /p choice=Enter [1/2]?: 
echo.
::echo choice = %choice%
echo.
echo.
if "%choice%" == "1" (
  echo entered 1: accuracy test
  build_accuracy
  echo.
  echo ********************************************
  echo.
  echo RUN EXECUTABLE FOR ACCURACY TEST
  echo     ex4_acc 1e-6 4 2
  ex4_acc 1e-6 4 2
  echo.
  echo.
  echo ********************************************
  echo.
  echo    TO REMOVE EXECUTABLE AND OBJECT FILES,
  echo    ENTER:
  echo         clean
  echo.
  echo    OTHERWISE RUN AS:
  echo         ex4_acc tol thrds1
  echo    OR
  echo         ex4_acc tol thrds1 thrds2
  echo.
  echo ********************************************
  echo.

) else (

  echo entered 2: time test
  build_times
  echo.
  echo ********************************************
  echo.
  echo RUN EXECUTABLE FOR TIME TEST
  echo     ex4_time 1e-12 1 5 5
  ex4_time 1e-12 1 5 5
  echo.
  echo.
  echo ********************************************
  echo.
  echo    TO REMOVE EXECUTABLE AND OBJECT FILES,
  echo    ENTER:
  echo         clean
  echo.
  echo    OTHERWISE RUN AS:
  echo         ex4_time tol jFUN NTval NXval
  echo.
  echo ********************************************
  echo.
)
