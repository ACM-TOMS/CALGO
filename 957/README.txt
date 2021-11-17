Since this software was constructed and tested it has been found
that the Matlab eigenvalue solver available since at least 2012
returns inaccurate eigenvectors which cause inaccuracies in
the double precision computed values of the repeated integral
of the coerror function.

This is still a problem at Matlab release 2016a. The directory
TestMatlab contains the sources that illustrate this problem;
the error in the last column of the file testRes.txt shows
the problem. This example has been taken from the book
     Walter Gautschi, "Orthogonal Polynomial in MATLAB:
     Exercises and Solutions", Software, Environments,
     Tools, SIAM, Philadelphia, 2016
and Demo 1.17 appears on pp.31--33. 
[The software is available from:
   http://www.siam.org/books/se26/
click on the "MATLAB files" link directly below the picture of
the book to download all the files associated with the book
including the files in TestMatlab]

********************************************************** 
This release has NOT been tested with a recent version of
Matlab; please report any problems (or, hopefully, fixes) to
t.r.hopkins@kent.ac.uk
********************************************************** 

Tim Hopkins
Algorithms Editor TOMS
June 14 2016
