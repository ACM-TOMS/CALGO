clear;

A1 = eye(2);
A2 = 2*eye(2);
A3 = 3*eye(2);
A4 = 4*eye(2);

A{1} = {[1 0 0 0],A1};
A{2} = {[1 0 0 0],A2};
A{3} = {[0 1 0 0],A3};
A{4} = {[0 0 0 1],A4};
poly = rolmipvar(A,'A',4,1);

dotpoly = diff(poly,'dpoly',[-1 2; -3 4; -8 6; -5 9])

