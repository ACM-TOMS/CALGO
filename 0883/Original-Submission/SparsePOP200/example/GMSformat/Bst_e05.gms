*  NLP written by GAMS Convert at 08/29/02 12:49:53
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         4       4       0       0       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         6       6       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        11       7       4       0
*
*  Solve m using NLP minimizing objvar;

Variables  x1,x2,x3,x4,x5,objvar;

Positive Variables x1,x2,x3,x4,x5;

Equations  e1,e2,e3,e4;


* e1.. 100000*x4 - 120*x1*(300 - x4) =E= 10000000;
e1.. 100000*x4 - 36000*x1 + 120*x1*x4 =E= 10000000;

* e2.. 100000*x5 - 80*x2*(400 - x5) - 100000*x4 =E= 0;
e2.. 100000*x5 - 32000*x2 + 80*x2*x5 - 100000*x4 =E= 0;

e3..  - 4000*x3 - 100000*x5 =E= -50000000;

e4..  - x1 - x2 - x3 + objvar =E= 0;

* set non default bounds

x1.up = 15834; 
x2.up = 36250; 
x3.up = 10000; 
x4.lo = 100; 
x4.up = 300; 
x5.lo = 100; 
x5.up = 400; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
