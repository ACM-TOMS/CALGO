*  NLP written by GAMS Convert at 09/01/02 13:58:14
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*        12       4       0       8       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         6       6       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        28      26       2       0
*
*  Solve m using NLP minimizing objvar;

***********************************************************
*
* To solve this problem by SparsePOP, set
*       param.relaxOrder = 2;
*
**********************************************************


Variables  x1,x2,objvar,x4,x5,x6;

Positive Variables x1,x2;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12;


e1..    x1 - 2*x2 =L= 100;

e2..  - 3*x1 - 4*x2 =L= -12;

e3..  - x1 - x2 =L= 100;

e4..  - x1 + 4*x2 =L= 100;

e5..  - x1 + 2*x2 =L= 18;

e6..    3*x1 + 4*x2 =L= 100;

e7..    x1 + x2 =L= 13;

e8..    x1 - 4*x2 =L= 8;

e9..    x1 - x4 =E= 0;

e10..    x1 - x2 - x5 =E= -10;

e11..    x1 + x2 - x6 =E= 6;

e12..  - x5*x6 + objvar - x4 =E= 0;

* set non default bounds

x1.up = 13; 
x2.up = 13; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
