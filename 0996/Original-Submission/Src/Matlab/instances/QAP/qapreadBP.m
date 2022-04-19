function [objPoly, I01, Icomp, relaxOrder, params] = qapreadBP(instance,lambda)
%% [objPoly, Icomp, I01, relaxOrder, params] = qapreadBP(filename,lambda)
%% Input
%%   instance : the instance name in the directory instances/QAP/qapdata
%%              for example, instance = 'chr12a'; 
%%   lambda   : the Lagrangian multiplier, for example lambda = 1.0e5
%% All instances in the directory instances/QAP/qapdata are from QAPLIB: 
%% P. Hahn and M. Anjos. QAPLIB ? a quadratic assignment problem library. 
%% http://www.seas.upenn.edu/qaplib.
%% For the DNN relaxation, see 
%% N. Ito, S. Kim, M. Kojima, A. Takeda, and K.-C. Toh.
%% A Sparse Doubly Nonnegative Relaxation of Polynomial Optimization Problems
%% with Binary, Box and Complementarity Constraints,
%% Research Rport B-48?, Department of Mathematical and Computing Sciences, 
%% Oh-Okayama, Meguro-ku, Tokyo 152-8552, March 2018. 

    filename = [instance '.dat'];
    
    fid = fopen(filename,'r');
    if (fid == -1); error('file cannot be opened'); end

    [datavec, ~] = fscanf(fid, '%c');
    fclose('all');

    linefeeds = strfind(datavec, char(10));
    datavec(linefeeds) = blanks(length(linefeeds));
    datavec = sscanf(datavec, '%f');

    nn = datavec(1);
    n2 = nn * nn;
    A = datavec(2:(n2 + 1));
    B = datavec((n2 + 2):(2 * n2 + 1));

    A = reshape(A, nn, nn);
    B = reshape(B, nn, nn);

    %% objPoly
    Q = kron(B, A);
    objPoly = quad2Poly(Q, [], []);

    %% complementarity
    % offdiag(X'X) = 0
    colCompPattern = repmat(speye(nn), nn, nn) - speye(n2);
    rowCompPattern = kron(speye(nn), ones(nn) - speye(nn));
    compPattern = colCompPattern | rowCompPattern;
    [rowIdx, colIdx] = find(triu(compPattern));
    compIdx = repmat(1:length(rowIdx), 2, 1);
    Icomp = sparse(compIdx, [rowIdx'; colIdx'], true, length(rowIdx), n2);

    %% binary
    % I01 = true(1, n2);
    I01 = false(1, n2);

    %% diag(X'X) = diag(I).
    C = [kron(ones(1, nn), speye(nn)); kron(speye(nn), ones(1, nn))];
    dd = ones(2 * nn, 1);
    penaltyPoly = quad2Poly(C' * C, -2 * (C' * dd), dd' * dd);

    %%%%%%%%%%
    penaltyPoly.coef = penaltyPoly.coef * lambda * norm(objPoly.coef) / norm(penaltyPoly.coef);
    objPoly = addPoly(objPoly, penaltyPoly);
    %%%%%%%%%%
    
    %%%%%%%%%%
    params.sparseSW = 0;
    relaxOrder = 1;
    params.delta = 0.1;
    params.delta1 = 1e-5;
    optVal = optvalQAP(instance); 
    params.UbdObjVal = optVal * 1.5;
    n = size(objPoly.supports, 2);
    params.UbdIX = 1 + sqrt(n);
    params.optVal = optVal; 
    %%%%%%%%%%
  
end