*  NLP written by GAMS Convert at 07/19/01 13:39:29
*  
*  Equation counts
*     Total       E       G       L       N       X
*        12       1       0      11       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        11      11       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*       112     105       7       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,objvar;

Positive Variables x1,x2,x3,x4,x5,x6,x7,x8,x9,x10;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12;


e1..  5*x1*x1 + 5*x2*x2 + 5*x3*x3 + 5*x4*x4 + 5*x5*x5 + 5*x6*x6
      + 5*x7*x7 + 20*x1 + 80*x2 + 20*x3 + 50*x4 + 60*x5 + 90*x6 - 10*x8
      - 10*x9 - 10*x10 + objvar =E= 0;

e2..  - 2*x1 - 6*x2 - x3 - 3*x5 - 3*x6 - 2*x7 - 6*x8 - 2*x9 - 2*x10 =L= -4;

e3..    6*x1 - 5*x2 + 8*x3 - 3*x4 + x6 + 3*x7 + 8*x8 + 9*x9 - 3*x10 =L= 22;

e4..  - 5*x1 + 6*x2 + 5*x3 + 3*x4 + 8*x5 - 8*x6 + 9*x7 + 2*x8 - 9*x10 =L= -6;

e5..    9*x1 + 5*x2 - 9*x4 + x5 - 8*x6 + 3*x7 - 9*x8 - 9*x9 - 3*x10 =L= -23;

e6..  - 8*x1 + 7*x2 - 4*x3 - 5*x4 - 9*x5 + x6 - 7*x7 - x8 + 3*x9 - 2*x10
      =L= -12;

e7..  - 7*x1 - 5*x2 - 2*x3 - 6*x5 - 6*x6 - 7*x7 - 6*x8 + 7*x9 + 7*x10 =L= -3;

e8..    x1 - 3*x2 - 3*x3 - 4*x4 - x5 - 4*x7 + x8 + 6*x9 =L= 1;

e9..    x1 - 2*x2 + 6*x3 + 9*x4 - 7*x6 + 9*x7 - 9*x8 - 6*x9 + 4*x10 =L= 12;

e10..  - 4*x1 + 6*x2 + 7*x3 + 2*x4 + 2*x5 + 6*x7 + 6*x8 - 7*x9 + 4*x10 =L= 15;

e11..    x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 =L= 9;

e12..  - x1 - x2 - x3 - x4 - x5 - x6 - x7 - x8 - x9 - x10 =L= -1;

* set non default bounds

x1.up = 1; 
x2.up = 1; 
x3.up = 1; 
x4.up = 1; 
x5.up = 1; 
x6.up = 1; 
x7.up = 1; 
x8.up = 1; 
x9.up = 1; 
x10.up = 1; 

* set non default levels

x1.l = 1; 
x2.l = 0.90755; 
x4.l = 1; 
x5.l = 0.71509; 
x6.l = 1; 
x8.l = 0.91698; 
x9.l = 1; 
x10.l = 1; 

* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
