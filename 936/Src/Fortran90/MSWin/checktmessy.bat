echo off
echo.
echo Compare contents of test output files newresult.d and result.d.
echo The latter is the original results from the algorithm distribution.
echo.
echo Values in this file from the algorithm distribution
echo should agree with those in newresult.d.
echo.
rem. MS-DOS fc = command for comparing file contents
fc ..\tmessy\result.d ..\tmessy\newresult.d /n /w
