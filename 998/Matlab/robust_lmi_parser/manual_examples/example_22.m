clear;

A{1} = {[2 0],eye(2)};
A{2} = {[1 1],3*eye(2)};
A{3} = {[0 2],5*eye(2)};

poly = rolmipvar(A,'A',2,2);
dpoly = partial(poly,'dAda')