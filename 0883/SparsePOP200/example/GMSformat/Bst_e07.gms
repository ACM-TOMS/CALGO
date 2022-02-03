*  NLP written by GAMS Convert at 08/29/02 12:49:53
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         8       6       0       2       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        11      11       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        32      25       7       0
*
*  Solve m using NLP minimizing objvar;


* Variables  x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,objvar;
* x10 = y10 -1; 
Variables  x1,x2,x3,x4,x5,x6,x7,x8,x9,y10,objvar;

Positive Variables x1,x2,x3,x4,x5,x6,x7,x8,x9,y10;

Equations  e1,e2,e3,e4,e5,e6,e7,e8;

e1..    x1 + x2 - x3 - x4 =E= 0;

e2..    x3 - x5 + x7 =E= 0;

e3..    x4 + x8 - x9 =E= 0;

e4..  - x6 + x7 + x8 =E= 0;

* e5.. x10*x3 - 2.5*x5 + 2*x7 =L= 0;
e5.. y10*x3 - x3 - 2.5*x5 + 2*x7 =L= 0;

* e6.. x10*x4 + 2*x8 - 1.5*x9 =L= 0;
e6.. y10*x4 - x4+ 2*x8 - 1.5*x9 =L= 0;

* e7..  - x3*x10 - x3*x4 + 3*x1 + x2 =E= 0;
e7..  - x3*y10 + x3 - x3*x4 + 3*x1 + x2 =E= 0;

e8..  - 6*x1 - 16*x2 + 9*x5 - 10*x6 + 15*x9 + objvar =E= 0;

* set non default bounds

x1.up = 300; 
x2.up = 300; 
x3.up = 100; 
x4.up = 200; 
x5.up = 100; 
x6.up = 300; 
x7.up = 100; 
x8.up = 200; 
x9.up = 200; 
* x10.lo = 1;
* x10.up = 3;
y10.up = 2;  

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
