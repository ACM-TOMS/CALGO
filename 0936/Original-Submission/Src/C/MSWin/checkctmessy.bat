echo off
echo.
echo Compare contents of test output files Results\result.d and result.d.
echo The latter is the original results from the algorithm distribution.
echo.
echo Values in this file from the algorithm distribution
echo should agree with those in Results\result.d
echo except for differences of -0. and +0.
echo.
echo fc = MS-DOS command for comparing file contents
echo.
fc ..\ctmessy\result.d ..\ctmessy\Results\result.d /n /w
