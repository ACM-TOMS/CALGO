RFEJER: Semi-Automatic Rational Fejer Quadrature

AUTHORS: Karl Deckers, Ahlem Mougaida, and Hedi Belhadjsalah

Software revision date: August 25, 2016

GENERAL DESCRIPTION
-------------------

RFEJER computes the nodes and weights in the n-point rational Fejer quadrature 
rule for a given sequence of n-1 real or complex conjugate poles outside the 
interval [-1,1]. In addition it also provides an approximation of the integral of 
an array-valued function FUN from -1 to 1, together with an estimation of the 
relative error, based on the computed nodes and weights. 'Semi' in 'semi-
automatic' refers to the fact that RFEJER does not use adaptive refinement of the 
interval, while the latter may be necessary to achieve an accuracy on the 
approximation of order machine precision. The algorithm, however, is implemented 
in such a way that it can easily be plugged into a larger program also containing 
an adaptive method.

===================================================================================

ALGORITHM PACKAGE
-----------------

The algorithm package contains:
* Main files: rfejer.m 
	with subroutines errW.m and transf.m
* Example files: exampleX.m (for X=1,...,7) 
	with subroutine nextpoint.m, output files exampleX.out (for X=1,...,7),
	and figures example7_PQ_figK.pdf (for PQ=11,12,31,32, and K=1,2,3)
* Test files: expaper.m, and expaperWG.m
* Third party supporting code: gauss.m, gauss_rational.m, lanczos.m, mcdis.m, 
	quadrat.m, rcheb.m, r_jacobi.m, r_mod.m, and stieltjes.m, 

All files need to be stored in the same directory. The third party supporting code
rcheb.m, which is used in the main file rfejer.m, has been updated to ensure 
compatibility with Octave (see further). The other third party supporting codes 
have not been updated and are only necessary for comparison with Gautschi's method 
in the example file example7.m and the test file expaperWG.m. 

===================================================================================

PORTABILITY
-------------

All files were initially developed in MATLAB R2011b, and have been tested 
afterwards to work properly in MATLAB R2013b, R2015a, R2015b, and R2016a as well. 
However, the file example7.m with second input argument 'intquad==1' does not work 
in MATLAB R2011b, since the MATLAB function integral.m was introduced later, in 
release 2012a. The output files exampleX.out, for X=1,...,7, were obtained by 
executing the corresponding MATLAB file exampleX.m in MATLAB R2011b on a notebook 
with a 2.26 GHz Intel Core 2 Duo P8400 processor. The figures 
example7_PQ_figK.pdf, for PQ=11,12,31,32, and K=1,2,3, were obtained by executing 
example7.m in MATLAB R2015b with an Intel(R) Core(TM) i7-3740QM CPU @2.70 GHz 
processor.

Compatibility with GNU Octave 4.0.3 requires:
1) the creation of a file named fcnchk.m with the following contents:
   _______________________
   function f=fcnchk(x, n)
   	f = x;
   end
   _______________________

2) removing '%' in 
   * example7.m, lines 134 and 136;
   * expaperWG.m, lines 116, 118, 140, and 142.

The third party supporting code rcheb.m has been updated to ensure compatibility 
with Octave. Whenever the initial version of rcheb.m is used instead, one should 
also:
3) remove '%' in rfejer.m, lines 198, 200, 231, and 233; and
4) make the following two small adaptations in rcheb.m:
   * line 353: 'bnd(2,i) == xo(i)' should be replaced by 'bnd(2,i) == xo(1,i)'
   * line 355: 'bnd(1,i) == xo(i)' should be replaced by 'bnd(1,i) == xo(1,i)'

Furthermore, when executing the test file expaper.m with input argument 'ex == 1, 
2, 3, or 8', messages of the form "warning: division by zero" can be avoided by 
using the alternatives from lines 325-353 instead on lines 38, 52, 76, and 110-113.

Finally, the file example7.m with second input argument 'intquad==1' does not 
work in Octave, since the latter has no built-in automatic integrator called
'integral'.

===================================================================================

DOCUMENTATION AND CONTACT INFORMATION
-------------------------------------

Deckers, K., Mougaida, A., and Belhadjsalah, H. 
Algorithm XXX: Extended Rational Fejer Quadrature Rules based on Chebyshev 
	       Orthogonal Rational Functions. 
ACM Trans. Math. Software (2016), 29 pages.

Corresponding author: Karl Deckers (karl.deckers@math.univ-lille1.fr)


