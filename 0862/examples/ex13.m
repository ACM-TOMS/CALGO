echo on
clear
% ----------------------------------------------------------------------
% Example 13: Tensor-tensor multiplication using tensor_as_matrix
% objects.
% ----------------------------------------------------------------------

% Create a 3 x 4 x 2 tensor
T = tensor(rand(3,4,2));

% Convert to a matrix
A = tensor_as_matrix(T,2)

% Compute the product of two matricized tensors
B = A' * A

% Or do the multiplication using ttt instead
C = ttt(T,T,2,2)


echo off
