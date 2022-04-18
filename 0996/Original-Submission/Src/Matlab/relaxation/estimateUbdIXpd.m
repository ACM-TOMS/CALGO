function UbdIX = estimateUbdIXpd(polyCone, Icomp, params)
%ESTIMATEUBDIXPD Estimate an uppre bound of <I,X> by a primal-dual algorithm
% Reference:
% N. Ito, S. Kim, M. Kojima, A. Takeda, and K.-C. Toh.
% A Sparse Doubly Nonnegative Relaxation of Polynomial Optimization Problems
% with Binary, Box and Complementarity Constraints,
% Research Rport B-48?, Department of Mathematical and Computing Sciences, 
% Oh-Okayama, Meguro-ku, Tokyo 152-8552, March 2018. 
%
% S. Iwata, and K. Nagano.
% Submodular Function Minimization under Covering Constraints
% 50th Annual IEEE Symposium on Foundations of Computer Science (FOCS), 2009


vS = polyCone.varStructure;
isDiag = vS.rowIdx == vS.colIdx;
isZero = polyCone.eq0Var(vS.momentic);
isNNZDiag = isDiag & ~isZero;
cardiB = sum(isNNZDiag);
NNZdiagic = vS.momentic(isNNZDiag);
diagSupports = vS.supports(NNZdiagic, :);

%[numV, numN] = size(diagSupports);
[numU, numN] = size(Icomp);
T = false(1, numN);
ST = any(Icomp(:, T), 2); %length==numU
z = zeros(1, numN);
y = zeros(numU, 1);
count = 0;
while ~all(ST)
    count = count + 1;
%     fprintf('PD update: %d', count);
    u = find(~ST, 1, 'first');
    Y = Icomp(u, :); %Y = Nu
    [alpha, X] = DinkelbachNewton(z, Y, diagSupports);
    y(u) = y(u) + alpha;
    z = z + alpha * Y;
    
    T = X;
    ST = any(Icomp(:, T), 2);
end
sumy = sum(y);
UB = full(sum(any(diagSupports(:, T), 2)));
UbdIX = floor(cardiB - sumy);

if params.printyesDNN >= 2
    fprintf('cardiB: %d, UbdIX:%d, sumy:%f, UB:%d \n', cardiB, UbdIX, sumy, UB);
end

end

function [alpha, Z] = DinkelbachNewton(z, Y, diagSupports)
    [m, n] = size(diagSupports);
    % generate capacity matrix for Dinic's algorithm
    numnode = 2+n+m;
    s = 1; t = numnode; N = (1:n) + 1; V = (1:m) + N(end);
    Capacity = spalloc(numnode, numnode, nnz(diagSupports) + m + n);
    Capacity(N, V) = spones(diagSupports)' * 1000000;
    Capacity(V, t) = 1;
    % init
    Z = true(1, n);
    val = -1;
    count = 0;
    while val < 0
        count = count + 1;
        alpha = ( sum(any(diagSupports(:, Z), 2)) - sum(z(Z)) ) / sum(Y & Z);
        Capacity(s, N) = (z + alpha * Y);
        [valc, S] = cutByDinic(Capacity, s, t);
        Z = ~S(N);
        val = valc - sum(z) - alpha * sum(Y); %adjust constant modular(ground set)
        alphapre = alpha;
    end
    %%
%     fprintf('    Newton update: %d\n', count);
end