echo on
clear
% ----------------------------------------------------------------------
% Example 16: Example of computing the series of n-mode products
% using tensor_as_matrix
% ----------------------------------------------------------------------

% Compute a random tensor T 
sizeT = [5,4,3,2];
T = tensor(rand(sizeT));

% Compute corresponding matrices
for n = 1:ndims(T)
    U{n} = rand(size(T,n),size(T,n));
end

% Use a tensor command to compute the
% tensor times the sequence of matrices
A = ttm(T,U);

% Manipulate the tensor instead in matrix form
% to compute the same result.
rdim = 2;
M = tensor_as_matrix(T, rdim, 'fc');
cdims = M.cindices;
B = U{rdim} * double(M);
B = B * kron(kron(U{cdims(3)},U{cdims(2)}),U{cdims(1)})';

% Reshape result to match A
B = reshape(B, [sizeT(rdim) sizeT(cdims)] );
[sdims sindx] = sort([rdim cdims]);
B = permute(B,sindx);

% Compare the results
dif = norm(A - tensor(B)) / norm(A)


echo off
