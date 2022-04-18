In order to test the package with the linear systems in the
directory 'data', the following instructions could be used: 

1) make

2) ./babdcr<data/example1a.txt


Available testing driver files are
example1a.txt
example1b.txt
example2a.txt
example2b.txt
example3a.txt
example3b.txt

Testing driver files with the name that ends with
'a' use the subroutines BABDCR_FACT and BABDCR_SOLV,
'b' use the subroutine BABDCR_FACTSOLV.



Files WRIGHTMATR, SWFIIIMATR and RANDMATR have boundary 
blocks stored after S_i and R_i blocks.


example1
data: WRIGHTMATR (coefficient matrix)
      WRIGHTRHS (right hand side) 
      WRIGHTSOL (solution) 
Wright example with NRWBLK=2 and NBLOKS=200. 
As shown in the paper [Wright 1993], WRIGHTMATR has the
coefficient matrix as expressed in equation (13), where 
NBLOKS=200 is the number of subintervals used in the 
discretization and NRWBLK=2 is the dimension of the 
differential system. Boundary blocks are placed in the 
last block row and the right hand side is set in order 
to obtain a solution with all the components equal to 1.


example2
data: SWFIIIMATR (coefficient matrix)
      SWFIIIRHS (right hand side)
      SWFIIISOL (solution)
Swirling flow III example with NRWBLK=6 and NBLOKS=512. 
This example represents the first linear system (in BABD 
form) arising from the solution of the nonlinear system
Phi(Y)=0 by means of the Newton method. The nonlinear
system describes the Mono Implicit Runge Kutta formula
applied to the Boundary Value Problem (see [Muir, Pancer, 
Jackson 2003])
y1'= y2
y2'= ( y1y4-y3y2 ) / epsilon
y3'= y4 
y4'= y5 
y5'= y6
y6'= -( y3y6+y1y2 ) / epsilon
-1<x<1
y3( -1 )= y3( 1 )= y4( -1 )= y4( 1 )= 0
y1( -1 )= -1, y1( 1 )= 1
This ODE is taken from the family of test problems 
Swirling Flow III (SWFIII) (see [Ascher, Mattheij, Russell 
1995]). The initial number of subintervals is set to 
NBLOKS=512 and the order of the differential system is 
NRWBLK=6. As initial guess Y^(0) of the Newton method (the 
r.h.s. of the considered system is -Phi(Y^(0))) we use
( Y^(0)_i )'= ( xi, 1, 0, 0, 0, 0 )
where xi= -1+ih,  i= 0,...,NBLOKS and h= 2/NBLOKS.


example3
data: RANDMATR (coefficient matrix)
      RANDRHS  (right hand side)
      RANDSOL  (solution)
random example with NRWBLK=3 and NBLOKS=10. 
