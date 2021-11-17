* example13.gms

* Variables x1, x2, x3,objvar;
Variables  x3,objvar;

* Positive Variables x1, x2;

Equations  e1,e2;

* minimize objvar = -2*x1 +3*x2 -2*x3
* minimize objvar = -2*1 +3*0 -2*x3
e1..    2*x3^2 + objvar =E= -2;

e2..    x3 =L= 3;

* x1 <= 2
* x1.up = 2;

* x2 <= 1
* x2.up = 1;

* end of example13.gms
