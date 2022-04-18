* example33.gms --- not square -- reduce --> square, feasible

Variables  x1,x2,x3,x4,x5,x6,objvar;

* Positive Variables x1,x2,x3,x4,x5;

Equations  e1,e2,e3,e4,e5,e6;

* minimize objvar = x1 + 2*x2 + 3*x3 +4*x4 +x5 +x6
e1..    x5 +x1 + 2*x2 + 3*x3 +4*x4 +x6 - objvar =E= 0;

e2..    x5 +x1 - x2 + x3 -x4 +x6 =E= 2; 

e3..    x1 + x2  =E= 1;

e4..    x2 + x3  =E= 2;

e5..    x3 + x4  =E= 3;

e6..    x4 + x5 +x6 =E= 4;

* e7..    x1 - x2 + x3 -x4 +x5 +x6 =E= 2; 

* Minimum solution: 
* 

* x1 <= 2
*x1.up = 2;

* x2 <= 1
* x2.up = 1;

* end of example33.gms
