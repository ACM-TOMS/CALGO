README FILE FOR COMPARISON WITH ALGORITHM AS 311

This directory (.../vauto/test/compare) contains programs used to compare
varma_llc.m with the 1997 program jam197.for of Jose A. Mauricio [1,2]. The
comparison was carried out by the program compare_with_as311.m, which calls
varma_llc.m and jam197.for (through a gateway function as311.c) for the same
test problems as used by test_varma_llc. All these three programs are included
here, as well as the output of the test which is in the file comparison.txt. To
repeat this analysis it is necessary to compile both the C-program and the
Fortran program and link them into a mex file. This could be carried out through
the Matlab command:
                         >> mex as311.c jam197.for

Under Windows it may be necessary to set up mex first (possibly with the aid of
Gnumex). If the compilation is successful the comparison would now be carried
out with:
                         >> compare_with_as311

which should produce results (more or less) identical to those contained in
comparison.txt.

NOTE: There is a small mistake in [2], log of 2*pi is only given with single
precision. This has been corrected in the version of jam197.for here.

[1] Mauricio, J. A. (1997). Algorithm AS 311: The exact likelihood function
    of a vector autoregressive moving average model. J. R. Statist. Soc. C:
    Applied Statistics, 46 no. 1, 157–171.

[2] Mauricio, J. A. (1997). JAM197.FOR (source code for [1]). Web page
    http://www.ucm.es/info/ecocuan/jam/jam197.for
