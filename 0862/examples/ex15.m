echo on
clear
% ----------------------------------------------------------------------
% Example 15: Two choices for converting a tensor to a matrix.
% ----------------------------------------------------------------------

% Let T be a 3 x 4 x 2 tensor
T = tensor(rand(3,4,2))

% Backward cyclic
A1 = tensor_as_matrix(T,2,'bc')

% Forward cyclic
A2 = tensor_as_matrix(T,2,'fc')

echo off