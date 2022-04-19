clear;

A1 = eye(2);
A2 = 2*eye(2);
B1 = 3*eye(2);
B2 = 4*eye(2);

A{1} = {[1 0],A1};
A{2} = {[0 1],A2};
B{1} = {[1 0],B1};
B{2} = {[0 1],B2};

polyA = rolmipvar(A,'A',2,1);
polyB = rolmipvar(B,'B',2,1);
polyB = fork(polyB,'B',1,2);

polyT = [polyA polyB; 
    polyB zeros(2,2)];