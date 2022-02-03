function B = hooi(A,R,kmax)
%HOOI Higher-order orthogonal iteration
%
%   B = HOOI(A,R) computes the best rank(R1,R2,..,Rn) approximation of
%   tensor A, according to the specified dimensions in vector R.
%   The result returned in B is stored as a TUCKER_TENSOR.
%

A = tensor(A);
N = ndims(A);

% Default value
if ~exist('kmax','var')
    kmax = 5;
end

% Compute an orthonormal basis for the dominant
% Rn-dimensional left singular subspace of
% A_(n) (1 <= n <= N). We store its transpose.
for n = 1:N
    [Un, S, V] = ...
	svds(double(tensor_as_matrix(A,n)), R(n));
    Ut{n} = Un';
end

% Iterate until convergence
for k = 1:kmax
    for n = 1:N
	Utilde = ttm(A, Ut, -n);
	% Maximize norm(Utilde x_n W') wrt W and
	% keeping orthonormality of W
	[Un,S,V] = ...
	    svds(double(tensor_as_matrix(Utilde, n)), R(n));
	Ut{n} = Un';
    end
end

% Create the core array
lambda = ttm(A, Ut);

% Create cell array containing U from the cell
% array containing its transpose
for n = 1:N
    U{n} = Ut{n}';
end

% Assemble the resulting tensor
B = tucker_tensor(lambda, U);
