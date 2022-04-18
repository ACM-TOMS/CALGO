*  NLP written by GAMS Convert at 08/31/02 19:03:58
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         8       1       7       0       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         9       9       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        23      15       8       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,x7,x8,objvar;

Equations  e1,e2,e3,e4,e5,e6,e7,e8;


e1..  - x1 + x2 =G= -1;

e2..  - x2 + x3 =G= -1.05;

e3..  - x3 + x4 =G= -1.1;

e4..  - x4 + x5 =G= -1.15;

e5..  - x5 + x6 =G= -1.2;

e6..  - x6 + x7 =G= -1.25;

e7..  - x7 + x8 =G= -1.3;

e8..  1.69*x1*x1 + 7*x1 + x1*x2 + 6*x2 + 2*x1*x3 + 5*x3 + 3*x1*x4 + 4*x4 + 4
     *x1*x5 + 3*x5 + 5*x1*x6 + 2*x6 + 6*x1*x7 + x7 + 7*x1*x8 + x2*x1 + 1.69*x2*
     x2 + x2*x3 + 2*x2*x4 + 3*x2*x5 + 4*x2*x6 + 5*x2*x7 + 6*x2*x8 + 2*x3*x1 + 
     x3*x2 + 1.69*x3*x3 + x3*x4 + 2*x3*x5 + 3*x3*x6 + 4*x3*x7 + 5*x3*x8 + 3*x4*
     x1 + 2*x4*x2 + x4*x3 + 1.69*x4*x4 + x4*x5 + 2*x4*x6 + 3*x4*x7 + 4*x4*x8 + 
     4*x5*x1 + 3*x5*x2 + 2*x5*x3 + x5*x4 + 1.69*x5*x5 + x5*x6 + 2*x5*x7 + 3*x5*
     x8 + 5*x6*x1 + 4*x6*x2 + 3*x6*x3 + 2*x6*x4 + x6*x5 + 1.69*x6*x6 + x6*x7 + 
     2*x6*x8 + 6*x7*x1 + 5*x7*x2 + 4*x7*x3 + 3*x7*x4 + 2*x7*x5 + x7*x6 + 1.69*
     x7*x7 + x7*x8 + 7*x8*x1 + 6*x8*x2 + 5*x8*x3 + 4*x8*x4 + 3*x8*x5 + 2*x8*x6
      + x8*x7 + 1.69*x8*x8 - objvar =E= 0;

* set non default bounds

x1.lo = -1; 
x1.up = 1; 
x2.lo = -2.1; 
x2.up = 2; 
x3.lo = -3.2; 
x3.up = 3; 
x4.lo = -4.3; 
x4.up = 4; 
x5.lo = -5.4; 
x5.up = 5; 
x6.lo = -6.5; 
x6.up = 6; 
x7.lo = -7.6; 
x7.up = 7; 
x8.lo = -8.7; 
x8.up = 8; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
