*  NLP written by GAMS Convert at 09/01/02 13:47:36
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         8       4       1       3       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         6       6       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        20      18       2       0
*
*  Solve m using NLP minimizing objvar;

***********************************************************
*
* To solve this problem by SparsePOP, set
*       param.relaxOrder = 2;
*
**********************************************************


Variables  x1,x2,x3,x4,x5,objvar;

Positive Variables x1;

Equations  e1,e2,e3,e4,e5,e6,e7,e8;


e1..    2*x1 + 3*x2 =G= 9;

e2..    3*x1 - x2 =L= 8;

e3..  - x1 + 2*x2 =L= 8;

e4..    x1 + 2*x2 =L= 12;

e5..  - x4*x5 - x3 + objvar =E= 0;

e6..    x1 - x3 =E= 0;

e7..    x1 - x2 - x4 =E= -5;

e8..    x1 + x2 - x5 =E= 1;

* set non default bounds

x1.up = 12; 
x2.lo = 3; 
x2.up = 6; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
