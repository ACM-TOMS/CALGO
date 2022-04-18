*  NLP written by GAMS Convert at 07/19/01 13:39:32
*  
*  Equation counts
*     Total       E       G       L       N       X
*         7       1       0       6       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         6       6       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        30       1      29       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,objvar;

Positive Variables  x1,x2,x3,x4,x5;

Equations  e1,e2,e3,e4,e5,e6,e7;


e1..  0.8356891*x1*x5 + 37.293239*x1 + 5.3578547*x3*x3 - objvar =E= 40792.141;

e2.. 0.0056858*x2*x5 - 0.0022053*x3*x5 + 0.0006262*x1*x4 =L= 6.665593;

e3.. 0.0022053*x3*x5 - 0.0056858*x2*x5 - 0.0006262*x1*x4 =L= 85.334407;

e4.. 0.0071317*x2*x5 + 0.0021813*x3*x3 + 0.0029955*x1*x2 =L= 29.48751;

e5.. -0.0071317*x2*x5 - 0.0021813*x3*x3 - 0.0029955*x1*x2 =L= -9.48751;

e6.. 0.0047026*x3*x5 + 0.0019085*x3*x4 + 0.0012547*x1*x3 =L= 15.599039;

e7.. -0.0047026*x3*x5 - 0.0019085*x3*x4 - 0.0012547*x1*x3 =L= -10.699039;

* set non default bounds

x1.lo = 78; 
x1.up = 102; 
x2.lo = 33; 
x2.up = 45; 
x3.lo = 27; 
x3.up = 45; 
x4.lo = 27; 
x4.up = 45; 
x5.lo = 27; 
x5.up = 45; 

* set non default levels

x3.l = 29.9953; 
x4.l = 45; 
x5.l = 36.7758; 

* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
