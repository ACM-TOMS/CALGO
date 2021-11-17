*  NLP written by GAMS Convert at 09/01/02 13:49:33
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         9       3       0       6       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         5       5       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        21      19       2       0
*
*  Solve m using NLP minimizing objvar;

***********************************************************
*
* To solve this problem by SparsePOP, set
*       param.relaxOrder = 2;
*
**********************************************************


Variables  x1,x2,objvar,x4,x5;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9;


e1..  - 4*x1 + x2 =L= 12;

e2..  - 4*x1 - 2*x2 =L= 12;

e3..    x1 - 2*x2 =L= 6;

e4..    x1 - x2 =L= 3;

e5..    x1 + x2 =L= 2;

e6..    2*x1 + x2 =L= 2;

e7..  - x4*x5 + objvar =E= 0;

e8..    x1 + x2 - x4 =E= 0;

e9..    x1 - x2 - x5 =E= 0;

* set non default bounds

x1.lo = -3; 
x1.up = 1; 
x2.lo = -3; 
x2.up = 4; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
