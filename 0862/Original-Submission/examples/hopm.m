function B = hopm(A,kmax)
%HOPM Higher-order power method
%
%   B = HOPM(A) computes the best rank-1 approximation of
%   tensor A using the power method.  The result returned in B is
%   stored as a CP_TENSOR.
%

A = tensor(A);
N = ndims(A);

% Default value
if ~exist('kmax','var')
    kmax = 5;
end

% Compute the dominant left singluar vectors
% of A_(n) (2 <= n <= N)
for n = 2:N
    [u{n}, lambda(n), V] = ...
        svds(double(tensor_as_matrix(A,n)), 1);
end

% Iterate until convergence
for k = 1:kmax
    for n = 1:N
        u{n} = ttv(A, u, -n);
        lambda(n) = norm(u{n});
        u{n} = double(u{n}./lambda(n));
    end
end

% Assemble the resulting tensor
B = cp_tensor(lambda(N), u);
