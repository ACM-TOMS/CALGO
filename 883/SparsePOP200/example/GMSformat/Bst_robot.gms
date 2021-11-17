*  NLP written by GAMS Convert at 09/03/02 17:54:06
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         9       9       0       0       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         9       9       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        25       9      16       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,x7,x8,objvar;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9;


e1.. 0.004731*x1*x3 - 0.1238*x1 - 0.3578*x2*x3 - 0.001637*x2 - 0.9338*x4 + x7
      =E= 0.3571;

e2.. 0.2238*x1*x3 + 0.2638*x1 + 0.7623*x2*x3 - 0.07745*x2 - 0.6734*x4 - x7
      =E= 0.6022;

e3.. x6*x8 + 0.3578*x1 + 0.004731*x2 =E= 0;

e4..  - 0.7623*x1 + 0.2238*x2 =E= -0.3461;

*e5.. POWER(x1,2) + POWER(x2,2) =E= 1;
e5.. x1*x1 + x2*x2 =E= 1;

*e6.. POWER(x3,2) + POWER(x4,2) =E= 1;
e6.. x3*x3 + x4*x4 =E= 1;

*e7.. POWER(x5,2) + POWER(x6,2) =E= 1;
e7.. x5*x5 + x6*x6 =E= 1;

*e8.. POWER(x7,2) + POWER(x8,2) =E= 1;
e8.. x7*x7 + x8*x8 =E= 1;

e9..    0 + objvar =E= 0;

* set non default bounds

x1.lo = -1; 
x1.up = 1; 
x2.lo = -1; 
x2.up = 1; 
x3.lo = -1; 
x3.up = 1; 
x4.lo = -1; 
x4.up = 1; 
x5.lo = -1; 
x5.up = 1; 
x6.lo = -1; 
x6.up = 1; 
x7.lo = -1; 
x7.up = 1; 
x8.lo = -1; 
x8.up = 1; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;