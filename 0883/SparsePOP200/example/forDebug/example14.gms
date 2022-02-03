* example14.gms --- fixed all variables, feasible.

Variables  x1,x2,x3,x4,x5,objvar;

Positive Variables x1, x2;

Equations  e1,e2,e3,e4,e5,e6;

* minimize objvar = -2*x1 +3*x2 -2*x3
e1..    2*x1 - 3*x2 + 2*x3 + objvar =E= 0;

e2..    x1 =E= 1;

e3..    x1 + x2 =E= 1;

e4..    x1 + x2 + x3 =E= 1;

e5..    x1 + x2 + x3 + x4 =E= 1;

e6..    x1 + x2 + x3 + x4 + x5 =E= 1;

* x1 <= 2
x1.up = 2;

* x2 <= 1
x2.up = 1;

* end of example10.gms
