*  NLP written by GAMS Convert at 08/31/02 19:06:06
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*        14       1       2      11       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        11      11       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*       116     106      10       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,objvar;

Positive Variables x1,x2,x3,x4,x5,x6,x7,x8,x9,x10;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14;


e1..    x1 + 7*x2 + 5*x3 + 5*x4 - 6*x6 - 3*x7 - 3*x8 + 5*x9 - 7*x10 =L= 80;

e2..  - 3*x1 + 3*x2 + 8*x3 + 7*x4 - 9*x5 - 7*x6 - 9*x7 + 8*x9 - 7*x10 =L= 57;

e3..    x1 + x3 + 3*x4 + 8*x5 + 9*x6 + 9*x8 - 7*x9 - 8*x10 =L= 92;

e4..  - x1 - 2*x2 + 2*x3 + 9*x5 + 5*x6 - 3*x7 + x8 - x9 - 5*x10 =L= 55;

e5..  - 5*x1 + 8*x2 - 8*x3 + 3*x5 + 4*x7 - 5*x8 - 2*x9 + 9*x10 =L= 76;

e6..    4*x1 - x2 + 6*x3 - 4*x4 - 7*x5 - 8*x6 - 7*x7 + 6*x8 - 2*x9 - 9*x10
      =L= 14;

e7..    7*x2 + 4*x3 + 9*x5 - 6*x8 - 5*x9 - 5*x10 =L= 47;

e8..  - 5*x1 - x2 + 7*x4 - x5 + 2*x6 + 5*x7 - 8*x8 - 5*x9 + 2*x10 =L= 51;

e9..  - 4*x1 - 7*x2 - 9*x4 + 2*x5 + 6*x6 - 9*x7 + x8 - 5*x9 =L= 36;

e10..  - 2*x1 + 6*x2 + 8*x4 - 6*x5 + 8*x6 + 8*x7 + 5*x8 + 2*x9 - 7*x10 =L= 92;

e11..    x1 + x2 + x3 - 2*x4 + x5 + x6 + x7 + 4*x8 + x9 + 3*x10 =L= 200;

e12..    x1 + x2 + x3 + x4 + x5 =G= 1;

e13..    x6 + x7 + x8 + x9 + x10 =G= 2;

e14..  x1*x6 - x1 - x6 + x2*x7 - 2*x2 - 2*x7 + x3*x8 - 3*x3 - 3*x8 + x4*x9
       - 4*x4 - 4*x9 + x5*x10 - 5*x5 - 5*x10 - objvar =E= 0;

* set non default bounds

x1.up = 100; 
x2.up = 100; 
x3.up = 100; 
x4.up = 100; 
x5.up = 100; 
x6.up = 100; 
x7.up = 100; 
x8.up = 100; 
x9.up = 100; 
x10.up = 100; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
