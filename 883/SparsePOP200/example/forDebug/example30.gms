* example30.gms

Variables  x1,x2,x3,x4,x5,objvar;

* Positive Variables x1,x2,x3,x4,x5;

Equations  e1,e2,e3,e4,e5,e6;

* minimize objvar = x1 + 2*x2 + 3*x3 +4*x4
e1..    x1 + 2*x2 + 3*x3 +4*x4 +5*x5 - objvar =E= 0;

e2..    x1 + x2  =E= 1;

e3..    x2 + x3  =E= 2;

e4..    x3 + x4  =E= 3;

e5..    x4 + x5  =E= 4;

e6..    x1 - x2 + x3 -x4 +x5 =E= 2; 

* Minimum solution: 
* 

* x1 <= 2
*x1.up = 2;

* x2 <= 1
* x2.up = 1;

* end of example30.gms
