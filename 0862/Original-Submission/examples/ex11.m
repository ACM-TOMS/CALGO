echo on
clear
% ----------------------------------------------------------------------
% Example 11: Computing the inner product with ttt or *
% ----------------------------------------------------------------------

% Create a 3 x 4 x 2 tensor
A = tensor(rand(3,4,2))

% Create a 3 x 4 x 2 tensor
B = tensor(rand(3,4,2))

% Tensor inner product of same-sized tensors
C = ttt(A,B,[1:3])

% Can also be computed using '*' notation
D = A * B;

% Results are equal
norm(C - D)

echo off