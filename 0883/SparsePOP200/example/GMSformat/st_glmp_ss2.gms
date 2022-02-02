*  NLP written by GAMS Convert at 09/01/02 14:01:00
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         9       4       1       4       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         6       6       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        22      20       2       0
*
*  Solve m using NLP minimizing objvar;

***********************************************************
*
* To solve this problem by SparsePOP, set
*       param.relaxOrder = 2;
*
**********************************************************


Variables  x1,x2,x3,x4,x5,objvar;

Positive Variables x1,x2;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9;


e1..    x1 - 2*x2 =L= 100;

e2..  - x1 - 2*x2 =L= 100;

e3..  - x1 + 2*x2 =G= 5;

e4..  - x1 + 2*x2 =L= 8;

e5..    x1 + 2*x2 =L= 12;

e6..    x1 - x3 =E= 0;

e7..    2*x1 - 3*x2 - x4 =E= -13;

e8..    x1 + x2 - x5 =E= 1;

e9..  - x4*x5 - x3 + objvar =E= 0;

* set non default bounds

x1.up = 7; 
x2.up = 6; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
