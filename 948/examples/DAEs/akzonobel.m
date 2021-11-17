function f = akzonobel(t,y)
% f = akzonobel(t,y) evaluates the chemical akzo nobel problem, used in 
% Section 5.3 of R. McKenzie, J. Pryce, N. Nedialkov, G. Tan 
% "DAESA User Guide".
%
% Copyright 2014 N. Nedialkov, J. Pryce, G. Tan

k1 = 18.7; k2 = 0.58;   k3 = 0.09; k4 = 0.42;
K  = 34.4; klA = 3.3; CO2 = 0.9;    H = 737;
Ks = 115.83;

r1  = k1*y(1)^4*sqrt(y(2));
r2  = k2*y(3)*y(4);
r3  = k2/K*y(1)*y(5);
r4  = k3*y(1)*y(4)^2;
r5  = k4*y(6)^2*sqrt(y(2));
Fin = klA*(CO2/H - y(2));

f(1) = -Dif(y(1),1) - 2.0*r1 + r2 - r3 - r4;
f(2) = -Dif(y(2),1) - 0.5*r1 - r4 - 0.5*r5 + Fin;
f(3) = -Dif(y(3),1) + r1 - r2 + r3;
f(4) = -Dif(y(4),1) - r2 + r3 - 2.0*r4;
f(5) = -Dif(y(5),1) + r2 - r3 + r5;
f(6) = Ks*y(1)*y(4) - y(6);

end
