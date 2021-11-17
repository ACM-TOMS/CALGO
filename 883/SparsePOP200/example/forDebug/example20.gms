* example20.gms --- feasible, optimal

Variables  x1,x2,x3,objvar;

Positive Variables x3;

Equations  e1,e2,e3, e4;

* minimize objvar = x1 + 2*x2 + x3
e1..    x1 + 2*x2 + x3 - objvar =E= 0;

e2..    x1 + 2*x2 =E= 1;

e3..    x1 + 2*x2 =E= 1;

e4..    x1 + 2*x2 =E= 1;

* x1 <= 2
*x1.up = 2;

* x2 <= 1
* x2.up = 1;

* end of example20.gms
