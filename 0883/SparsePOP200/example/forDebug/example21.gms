* example21.gms --- feasible, optimal

Variables  x1,x2,x3,x4,x5,objvar;

* Positive Variables x1,x2,x3,x4,x5;

Equations  e1,e2,e3, e4, e5;
* , e4;

* minimize objvar = x1 + 2*x2 + x3
* e1..    x1 + 6*x2 + 3*x3 + x4^2 +x5^2 - objvar =E= 0;
e1..    x1 + 2*x4 + 3*x3 + x2^2 +x5^2 - objvar =E= 0;

e2..    x1 + 2*x4 + 3*x3 -x2 =E= 1;

e3..    x1 + 2*x4 + 3*x3 -x5 =E= 4;

e4..    x1 + 2*x4 + 3*x3 -x5 =E= 4;

e5..    2*x1 + 4*x4 + 6*x3 -2*x2 =E= 2;

* Minimum solution: 
* 

* x1 <= 2
*x1.up = 2;

* x2 <= 1
* x2.up = 1;

* end of example21.gms
