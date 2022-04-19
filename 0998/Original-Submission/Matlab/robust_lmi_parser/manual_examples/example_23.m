clear;

A{1} = eye(2);
A{2} = 3*eye(2);

polyA = rolmipvar(A,'A',2,1);
polyout = discshift(polyA,2,[-0.4 0.4; -0.6 0.6]);

polyout{1}
polyout{2}
polyout{3}