function lambda=DirichletEig(S, M)
% compute the first eigenvalue for -\Delta 
lambda=eigs(S, M, 1, 'SM'); 

