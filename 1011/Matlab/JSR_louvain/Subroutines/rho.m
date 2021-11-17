function rhoo = rho(M)
% ANS = RHO(M)
% 
% Computes the spectral radius of matrix M
% if M is a cell of matrices, ans is a row vector with 
% the spectral radius of each M{i}

opts.disp = 0;

if(iscell(M))
    m = length(M);
    rhoo = zeros(1,m);

    for i = 1:m
        if issparse(M{i})
            rhoo(i) = eigs(M{i},1,'LM',opts);
        else
            rhoo(i) = max(abs(eig(M{i})));
        end
    end
else
    if (issparse(M))
        rhoo = eigs(M,1,'LM',opts);
    else
        rhoo = max(abs(eig(M)));
    end
end
end