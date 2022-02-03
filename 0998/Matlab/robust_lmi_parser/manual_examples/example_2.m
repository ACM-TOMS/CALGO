clear;

A1 = eye(2);
A2 = 2*eye(2);
A3 = 3*eye(2);
A4 = 4*eye(2);
A5 = 5*eye(2);
A6 = 6*eye(2);
A7 = 7*eye(2);
A8 = 8*eye(2);
A9 = 9*eye(2);

A{1} = {[2 0],[1 0 0],A1};
A{2} = {[1 1],[1 0 0],A2};
A{3} = {[0 2],[1 0 0],A3};
A{4} = {[2 0],[0 1 0],A4};
A{5} = {[1 1],[0 1 0],A5};
A{6} = {[0 2],[0 1 0],A6};
A{7} = {[2 0],[0 0 1],A7};
A{8} = {[1 1],[0 0 1],A8};
A{9} = {[0 2],[0 0 1],A9};
poly_A = rolmipvar(A,'A',[2 3],[2 1]);

%Retrieve monomial A4
mon_A4 = poly_A([2 0],[0 1 0]);
mon_A4{1}


%Setting the monomial {[2 0],[0 1 0]} as a 3x3 identity matrix:
poly_A([2 0],[0 1 0]) = eye(3);
mon_A4 = poly_A([2 0],[0 1 0]);
mon_A4{1}