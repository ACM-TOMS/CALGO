clear;

A1 = eye(2);
A2 = 2*eye(2);
A3 = 3*eye(2);
A4 = 4*eye(2);

A{1} = {[1 0],[1 0],A1};
A{2} = {[0 1],[1 0],A2};
A{3} = {[1 0],[0 1],A3};
A{4} = {[0 1],[0 1],A4};
poly_A = rolmipvar(A,'A',[2 2],[1 1])

poly_Afork = fork(poly_A,'Afork')