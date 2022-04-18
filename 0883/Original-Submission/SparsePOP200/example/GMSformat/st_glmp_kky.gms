*  NLP written by GAMS Convert at 09/01/02 13:54:19
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*        14       6       0       8       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         8       8       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        37      33       4       0
*
*  Solve m using NLP minimizing objvar;

***********************************************************
*
* To solve this problem by SparsePOP, set
*       param.relaxOrder = 2;
*
**********************************************************


Variables  x1,x2,x3,x4,x5,x6,x7,objvar;

Positive Variables x1,x2;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14;


e1..  - 5*x1 + 8*x2 =L= 24;

e2..  - 5*x1 - 8*x2 =L= 100;

e3..  - 6*x1 + 3*x2 =L= 100;

e4..  - 4*x1 - 5*x2 =L= -10;

e5..    5*x1 - 8*x2 =L= 100;

e6..    5*x1 + 8*x2 =L= 44;

e7..    6*x1 - 3*x2 =L= 15;

e8..    4*x1 + 5*x2 =L= 100;

e9..  - (x4*x5 + x6*x7) - x3 + objvar =E= 0;

e10..    3*x1 - 4*x2 - x3 =E= 0;

e11..    x1 + 2*x2 - x4 =E= 1.5;

e12..    2*x1 - x2 - x5 =E= -4;

e13..    x1 - 2*x2 - x6 =E= -8.5;

e14..    2*x1 + x2 - x7 =E= 1;

* set non default bounds

x1.up = 10; 
x2.up = 10; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
