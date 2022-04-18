* exampleInfeaisble.gms

Variables  x1,x2,x3,objvar;

Positive Variables x3;

Equations  e1,e2,e3,e4;

* minimize objvar = x1 +3*x2 + 2*x3
e1..    x1 +3*x2 + 2*x3 - objvar =E= 0;

e2..    x1 + x2 =G= 7;

e3..    x2 + x3 =G= 9;

e4..    x1 + x3 =G= 3;

* x1 <= 2
x1.up = 2;

* x2 <= 1
x2.up = 1;

x3.up = 1; 

* end of example11.gms
