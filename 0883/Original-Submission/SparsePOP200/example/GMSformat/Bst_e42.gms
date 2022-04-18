*  NLP written by GAMS Convert at 08/29/02 12:49:54
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         3       3       0       0       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         8       8       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        11       7       4       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,x7,objvar;

Positive Variables x1,x2,x3,x4,x7;

Equations  e1,e2,e3;


e1.. 51.5712*x3*x5 + 20.7324*x3 - 25.7856*x5 + 14.9251*x3*x7 - 22.2988*x7 - 
     10.1409*x6*x7 + 42.3401*x6 - x1 + x2 - 6.6425*x4 =E= -72.82;

e2..    x3 =E= 1;

e3..  - x1 - x2 + objvar =E= 0;

* set non default bounds

x3.up = 1; 
x4.up = 1; 
x5.lo = -1; 
x5.up = 1; 
x6.lo = -1; 
x6.up = 1; 
x7.up = 3; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
