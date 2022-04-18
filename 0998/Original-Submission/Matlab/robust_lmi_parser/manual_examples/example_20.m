clear;

A1 = eye(2);
A2 = 2*eye(2);
A3 = 3*eye(2);
A4 = 4*eye(2);
A5 = 5*eye(2);
A6 = 6*eye(2);

A{1} = {[1 0 0],[1 0],A1};
A{2} = {[1 0 0],[0 1],A2};
A{3} = {[0 1 0],[1 0],A3};
A{4} = {[0 1 0],[0 1],A4};
A{5} = {[0 0 1],[1 0],A5};
A{6} = {[0 0 1],[0 1],A6};
poly = rolmipvar(A,'A',[3 2],[1 1]);

dotbounds{1} = [-1 2; -3 4; -8 6];
dotbounds{2} = [-4 1; -6 9];
dotpoly = diff(poly,'dpoly',dotbounds)

