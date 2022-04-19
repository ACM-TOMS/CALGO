*  NLP written by GAMS Convert at 08/29/02 12:53:30
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         6       5       0       1       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         7       7       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        17       7      10       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,objvar;

Positive Variables x1,x2,x3,x4;

Equations  e1,e2,e3,e4,e5,e6;


e1.. 0.09755988*x1*x5 + x1 =E= 1;

e2.. 0.0965842812*x2*x6 + x2 - x1 =E= 0;

e3.. 0.0391908*x3*x5 + x3 + x1 =E= 1;

e4.. 0.03527172*x4*x6 + x4 - x1 + x2 - x3 =E= 0;

e5.. x5**0.5 + x6**0.5 =L= 4;

e6..    x4 + objvar =E= 0;

* set non default bounds

x1.up = 1; 
x2.up = 1; 
x3.up = 1; 
x4.up = 1; 
x5.lo = 1E-5; 
x5.up = 16; 
x6.lo = 1E-5;
x6.up = 16; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
