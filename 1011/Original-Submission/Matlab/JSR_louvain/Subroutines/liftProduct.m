function P = liftProduct(M, k)

% LIFTPRODUCT Generation of the set of all k-products of matrices of a set
%
%    P = liftProduct(M, k) returns a cell array of all products of less or equal k matrices
%      taken from the set M. The same matrix can appear several times in the
%      product, and P may contain several identical products

n = length(M);
m = size(M{1}, 1);
P = liftProductRecursive(cell(n^k, 1), M, n, eye(m), k-1, 0);



function P = liftProductRecursive(P, M, n, prod, k, shift)

if (k <= 0),
    for i = 1:n,
        P{shift+i, 1} = prod*M{i};
    end
    return
end

for i = 1:n,
    P = liftProductRecursive(P, M, n, prod*M{i}, k-1, shift+(i-1)*n^k);
end