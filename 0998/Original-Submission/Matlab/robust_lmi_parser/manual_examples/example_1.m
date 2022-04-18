clear;

A1 = eye(2);
A2 = 2*eye(2);
A3 = 3*eye(2);

A = [A1 A2 A3];
poly_A = rolmipvar(A,'A',3,1);

%Retrieve monomial A2
mon_A2 = poly_A([0 1 0]);
mon_A2{1}