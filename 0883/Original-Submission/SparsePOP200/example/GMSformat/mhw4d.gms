*  NLP written by GAMS Convert at 07/26/01 11:56:47
*  
*  Equation counts
*     Total       E       G       L       N       X
*         4       4       0       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         6       6       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        14       4      10       0
*
*  Solve m using NLP minimizing objvar;


Variables  objvar,x2,x3,x4,x5,x6;

Equations  e1,e2,e3,e4;


e1..  - ((x2 - 1)*(x2 - 1) + (x2 - x3)*(x2 - x3) + (x3 - x4)*(x3 - x4)*(x3 - x4) + (x4 - x5)*(x4 - x5)*(x4 - x5)*(x4 - x5) + 
     (x5-x6)*(x5 - x6)*(x5-x6)*(x5 - x6)) + objvar =E= 0;

e2.. x3*x3 + x4*x4*x4 + x2 =E= 6.24264068711929;

e3..  - x4*x4 + x3 + x5 =E= 0.82842712474619;

e4.. x2*x6 =E= 2;

* set non default bounds


* set non default levels

x2.l = -1; 
x3.l = 2; 
x4.l = 1; 
x5.l = -2; 
x6.l = -2; 

* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
