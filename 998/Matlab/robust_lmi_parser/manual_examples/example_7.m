clear;

A0 = 0.1*eye(2);
A1 = eye(2);
A2 = 2*eye(2);
A3 = 3*eye(2);

A{1} = {[0 0],A0};
A{2} = {[1 0],A1};
A{3} = {[0 1],A2};
A{4} = {[2 1],A3};

poly_A = rolmipvar(A,'A',[-2 3; -4 8]);
