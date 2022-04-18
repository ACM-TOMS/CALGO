echo off
echo.
echo Compare contents of ...\newresults\err0.d, mes0.d and odd1.d with
echo those found in ...\Driver\Results of the algorithm distribution.
echo.
echo With 2 OpenMP threads, contents of files err0.d, mes0.d and odd1.d from
echo the algorithm distribution should agree with those in ...\newresults.
echo.
rem. MS-DOS fc = command for comparing file contents
fc err0.d newresults\err0.d /n /w
fc mes0.d newresults\mes0.d /n /w 
fc odd1.d newresults\odd1.d /n /w
