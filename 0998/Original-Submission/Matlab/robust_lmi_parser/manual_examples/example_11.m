clear;

A1 = eye(2);
A2 = 2*eye(2);
A = [A1 A2];
poly_A = rolmipvar(A,'A',2,1);

poly_Afork = fork(poly_A,'Afork')