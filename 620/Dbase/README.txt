
         REFERENCES AND KEYWORDS FOR ACM-CALGO ALGORITHMS

       Prepared Oct. 11, 1984 and modified Mar. 1, 1986 by
         J. R. Rice & R. J. Hanson

       Revised July 27, 1989 by Tim Hopkins and David Morse

       Revised June 9, 1992 and August 4, 1994 by Robert Renka

       Annual updates will be made without comment

This file contains compressed references for the set of algorithms
contained in Collected Algorithms From ACM.  Algorithms 1-492 were
published in Communications of the ACM.  With one exception, Algor-
ithms 493-701 appeared in ACM Transactions on Mathematical Software.
The exception is Algorithm 568, which was published in ACM Transac-
tions on Programming Languages and Systems.  Algorithm numbers 171
and 172 are not included because no algorithm with those numbers was
published.

The entry for each algorithm consists of either four or five records
depending on whether there have been any published remarks.  Lines
are restricted to 80 characters, and records requiring more than one
line include continuation lines which are distinguished by a + in
the first column.  Only the first record and continuation lines begin
in character column one.  Any mathematical notation used within the
algorithm title and any accents in an author's name are indicated as
defined in TEX.  Also, all letters in the title which must remain
capitalized in a printed version of the reference are enclosed in
braces.

The first record specifies the algorithm number, journal in which the
algorithm was published, beginning page number, ending page number (or
0 if there is only one page), volume number, issue number, month and
year of publication, modified SHARE classification, and language (F =
Fortran, F90 = Fortran 90, A60 = Algol 60, PAS = Pascal, PLI = PL1, M =
Matlab, L = Lisp, R = Ratfor, N = None) in which the algorithm was
implemented.  Some of the later algorithms have the associated gams
classification as the last field.  The second and third records contain
the authors' names and the title of the algorithm, respectively.  The
fourth record contains keywords separated by semicolons, and the fifth
record, if any, describes published remarks associated with the
algorithm.  Each remark, if any, is terminated by a semicolon, and
includes the following fields, separated by commas:  the type (Remark
or Certification), the journal in which the remark was published, the
page range (a single page number or a pair of numbers separated by --),
the volume number, the issue number, the month, the year, and the
author.  As an example, the following entry is for algorithm 487:

487    cacm  703  704 17 12  December 1974 s14   F
 J. Pomeranz;
 Exact Cumulative Distribution of the {K}olmogorov-{S}mirnov Statistic for
+ Small Samples
 goodness-of-fit testing;k-s statistic;k-s test;Kolmogorov-Smirnov test;
 R,toms,111,2,1,March,1976,J. Pomeranz;
+R,toms,285--294,3,3,September,1977,R. Kallman;

The first line should be interpreted as 'ACM CALGO Algorithm 487 appeared
in Comm. ACM, Volume 17, Number 12, December 1974, pages 703-704'.  The
algorithm was implemented in Fortran, and the SHARE classification is S14.
The title spans two lines and contains two characters which must remain
in upper case.  The second remark appeared as a Remark in ACM TOMS, Volume
3, Number 3, September 1977, pages 285-294.  The author is R. Kallman.
