function P = tliftproduct(M, k)
% P = tliftproduct(M, k) 
% Computes all products of length less or equal k.
% The same matrix can appear several times in the product. Also returns the empty product eye()
%
% Input:
%   M       cell array of matrices
%   k       integer, maximal length of products
%
% Output:
%   P       all products of length less or equal k
%
% E.g.: tliftproduct({2 3 5}, 3)
%
% Written for JSR-toolbox, Jungers et al,
% Modified by tommsch, 2018

n = length(M);
m = size(M{1}, 1);
P = cell(k+1,1);
dim = size(M{1},1);
P{1} = eye(dim);
for i = 1:k
    P{i+1} = tliftproduct_recursive(cell(n^i, 1), M, n, eye(m), i-1, 0); end;
P = vertcat(P{:}).';


end

function P = tliftproduct_recursive(P, M, n, prod, k, shift)

if (k <= 0)
    for i = 1:n
        P{shift+i, 1} = prod*M{i}; end;
    return; end;

for i = 1:n
    P = tliftproduct_recursive(P, M, n, prod*M{i}, k-1, shift+(i-1)*n^k); end;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 