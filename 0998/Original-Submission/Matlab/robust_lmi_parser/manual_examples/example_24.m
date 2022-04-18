clear;

A{1} = {[1 0],[1 0],eye(2)};
A{2} = {[1 0],[0 1],3*eye(2)};
A{3} = {[0 1],[1 0],5*eye(2)};
A{4} = {[0 1],[0 1],7*eye(2)};
polyA = rolmipvar(A,'A',[2 2],[1 1]);

%Generates A(\alpha(k+1),\beta(k))
polyaux = discshift(polyA,1,[-0.4 0.4; -0.6 0.6],1,3);

%Applies discshift on A(\alpha(k+1),\beta(k)) to generate
%A(\alpha(k+1),\beta(k+1))
polyout = discshift(polyaux{1},1,[-0.4 0.4; -0.6 0.6],2,4);
polyout{1}