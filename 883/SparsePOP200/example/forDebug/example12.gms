* example12.gms --- fixed 2 variables ---> minimization in x3, unbounded. 

Variables  x1,x2,x3,objvar;

Positive Variables x1, x2;

Equations  e1,e2,e3,e4;

* minimize objvar = -2*x1 +3*x2 -2*x3
e1..    2*x1 - 3*x2 + 2*x3^2 + objvar =E= 0;

e2..    x1 =E= 1;

e3..    x2 =E= 0;

e4..    x3 =L= 3;

* x1 <= 2
x1.up = 2;

* x2 <= 1
x2.up = 1;

* end of example12.gms
