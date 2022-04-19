*  NLP written by GAMS Convert at 07/19/01 13:39:30
*  
*  Equation counts
*     Total       E       G       L       N       X
*        11      11       0       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        25      25       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        73      49      24       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19
          ,x20,x21,x22,x23,x24,objvar;

Positive Variables x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17
          ,x18,x19,x20,x21,x22,x23,x24;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11;


e1..  300*x1 - 7*x1*x1 - 4*x2*x2 + 270*x2 - 6*x3*x3 + 460*x3 - 8*x4*x4 + 800
     *x4 - 12*x5*x5 + 740*x5 - 9*x6*x6 + 600*x6 - 14*x7*x7 + 540*x7 - 7*x8*x8
      + 380*x8 - 13*x9*x9 + 300*x9 - 12*x10*x10 + 490*x10 - 8*x11*x11 + 380*x11
      - 4*x12*x12 + 760*x12 - 7*x13*x13 + 430*x13 - 9*x14*x14 + 250*x14 - 16*
     x15*x15 + 390*x15 - 8*x16*x16 + 600*x16 - 4*x17*x17 + 210*x17 - 10*x18*x18
      + 830*x18 - 21*x19*x19 + 470*x19 - 13*x20*x20 + 680*x20 - 17*x21*x21 + 
     360*x21 - 9*x22*x22 + 290*x22 - 8*x23*x23 + 400*x23 - 4*x24*x24 + 310*x24
     - objvar =E= 0;

e2..    x1 + x2 + x3 + x4 =E= 8;

e3..    x5 + x6 + x7 + x8 =E= 24;

e4..    x9 + x10 + x11 + x12 =E= 20;

e5..    x13 + x14 + x15 + x16 =E= 24;

e6..    x17 + x18 + x19 + x20 =E= 16;

e7..    x21 + x22 + x23 + x24 =E= 12;

e8..    x1 + x5 + x9 + x13 + x17 + x21 =E= 29;

e9..    x2 + x6 + x10 + x14 + x18 + x22 =E= 41;

e10..    x3 + x7 + x11 + x15 + x19 + x23 =E= 13;

e11..    x4 + x8 + x12 + x16 + x20 + x24 =E= 21;

* set non default bounds

x1.up = 24; 
x2.up = 24; 
x3.up = 24; 
x4.up = 24; 
x5.up = 24; 
x6.up = 24; 
x7.up = 24; 
x8.up = 24; 
x9.up = 24; 
x10.up = 24; 
x11.up = 24; 
x12.up = 24; 
x13.up = 24; 
x14.up = 24; 
x15.up = 24; 
x16.up = 24; 
x17.up = 24; 
x18.up = 24; 
x19.up = 24; 
x20.up = 24; 
x21.up = 24; 
x22.up = 24; 
x23.up = 24; 
x24.up = 24; 

* set non default levels

x1.l = 6; 
x2.l = 2; 
x6.l = 3; 
x8.l = 21; 
x9.l = 20; 
x13.l = 24; 
x17.l = 3; 
x19.l = 13; 
x22.l = 12; 

* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
