*  NLP written by GAMS Convert at 07/19/01 13:40:18
*  
*  Equation counts
*     Total       E       G       L       N       X
*        10      10       0       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        11      11       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        26      18       8       0
*
*  Solve m using NLP minimizing objvar;


Variables  objvar,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11;

Positive Variables x2,x3,x4,x5,x6,x7,x8,x9,x10,x11;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10;


e1..  - objvar + x2 - 4*x3 =E= 0;
* e1.. x6 + objvar =E= 0; 

e2..  - 2*x2 + x3 + x4 =E= 0;

e3..    2*x2 + 5*x3 + x5 =E= 108;

e4..    2*x2 - 3*x3 + x6 =E= -4;

e5..  - x3 + x7 =E= 0;

e6.. x8*x4 =E= 0;

e7.. x9*x5 =E= 0;

e8.. x10*x6 =E= 0;
* e8.. x6 =E= 0; 

e9.. x11*x7 =E= 0;
* e9.. x11 =E= 0; 

e10..    x8 + 5*x9 - 3*x10 - x11 =E= -1;

* set non default bounds

*
* Lower and upper bounds of variables are updated by SparsePOP
*

x2.up = 19;
x3.up = 14; 
x4.up = 24; 
x5.up = 96; 
x6.up = 0.0016; 
x7.up = 140; 
x8.up = 200; 
x9.up = 120; 
x10.up = 200; 
x11.up = 0.0004; 

* x4.up = 200; 
* x5.up = 200; 
* x6.up = 200; 
* x7.up = 200; 
* x8.up = 200; 
* x9.up = 200; 
* x10.up = 200; 
* x11.up = 200; 

* set non default levels


* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;

